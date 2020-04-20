## IMPACT(COVID-19) was developed by Chris Kypridemos with contributions from Roberta
## Piroddi and Alex Alexiou
##
## Copyright (C) 2020 University of Liverpool, Chris Kypridemos
##
## IMPACT(COVID-19) is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.


# Simulation parameters -------------------------------------------------------
set.seed(42)
sim_horizon <- 50L # additional days beyond the simulation horizon of the transmission model
probabilities <- c(0.5, 0.025, 0.975, 0.1, 0.9, 0.2, 0.8) # uncertainty intervals for the outputs
n_cores <- parallel::detectCores()/2L # Number of available cores for parallelisation
show_plots <- FALSE
delete_previous_results <- TRUE
collect_indiv_traj <- FALSE # Writes to disk individual trajectories. Huge file!!!

# Load required packages
library(simmer)
library(simmer.bricks)
library(simmer.plot)
library(parallel)
library(Rcpp) # Make sure Rtools are installed in Windows
library(scales) # for percent formatting
library(ggplot2)

library(data.table)

sourceCpp("functions.cpp", cacheDir = normalizePath("./cache/"))
source("./prepare_inputs.R")
sim_horizon <- infected[, max(date) - min(date)] + 1L + sim_horizon

# Set resource utilisation ----------------------------------------------------

# Time (in days) from the onset of symptom until healthcare utilisation
time_to_seek_healthcare <- function(n) rgamma(n, 1.8, 0.43) # median 5 days from the onset of symptoms

# Assume 3 hours median stay in A&E
AE_utilisation <- function() rbinom(1L, 24L, 3/24)/24 # This parameter is not important

# Assume 6 days mean level 1 bed utilisation for those that they will never need ICU
bed1_utilisation <- function()  runif(1L, 4, 8) # 6 days mean stay, from 4 to 8 days

# Assume 1 days mean level 1 bed utilisation before escalation to ICU
bed_preICU_utilisation <- function() rpois(1L, 1) # 6 days mean stay, from 4 to 8 days

# Assume 12 days mean ICU bed utilisation before escalation to ICU
ICU_utilisation <- function() runif(1L, 10, 14) # 12 days mean stay, from 10 to 14 days

# Assume 0 days mean level 1 bed utilisation after de-escalation from ICU
bed_postICU_utilisation <- function() 0.00001 # occupies level 1 bed just for a few seconds



# Set trajectories ------------------------------------------------------------

patient <- trajectory("patients' path") %>%
  # Patient arrives at A&E and gets assessed
  visit("A&E", function()
    AE_utilisation()) %>%

  # Discharge to community or admit to hospital
  branch(
    function()
      get_traj_hospitalisation(get_attribute(hospital, "agegroup")),
    continue = c(TRUE), # Nowhere to continue to. TRUE for safety
    trajectory() %>% # discharge from A&E
      visit("A&E_discharge", sim_horizon + 1L),
    trajectory() %>% # Admit to hospital

      # Decide whether it will require ICU treatment
      branch(
        # Patient only occupies a level 1 bed or eventually admitted to ICU
        function()
          get_traj_icu(get_attribute(hospital, "agegroup")),
        continue = c(TRUE), # Nowhere to continue to

        # Option 1: patient occupies level 1 bed exclusively
        trajectory() %>%
          visit("bed1", function() bed1_utilisation()) %>%

          # decide whether patient dies in hospital or discharged to community
          branch(
            function()
              get_traj_nonicu_death(get_attribute(hospital, "agegroup")),
            continue = c(FALSE, FALSE),
            trajectory() %>% # discharge
              visit("bed1_discharge", sim_horizon + 1L),
            trajectory() %>% # in-hospital death
              visit("in-hospital_death", sim_horizon + 1L)
          ),

        # Option 2: patient occupies a level 1 bed then an ICU bed and then a level 1 bed
        # again. Patients can die in ICU or after de-escalation but not before
        # ICU as this would mean they had never been admitted to ICU which is a
        # requirement for this pathway
        trajectory() %>% #
          visit("bed1", function()
            bed_preICU_utilisation()) %>%
          visit("ICU", function()
            ICU_utilisation()) %>%

          # decide whether patient dies in ICU
          branch(
            function()
              get_traj_icu_death(get_attribute(hospital, "agegroup")),
            continue = c(FALSE, TRUE),
            trajectory() %>% # death in ICU
              visit("in-hospital_death", sim_horizon + 1L),
            trajectory() %>% # Survives ICU and moves to level 1 bed
              visit("bed1", function()
                bed_postICU_utilisation()) %>%

              # decide whether patient dies post-ICU
              branch(
                function()
                  get_traj_nonicu_death(get_attribute(hospital, "agegroup")),
                continue = c(FALSE, FALSE),
                trajectory() %>% # discharge
                  visit("ICU_discharge", sim_horizon + 1L),
                trajectory() %>% # in-hospital death post ICU
                  visit("in-hospital_death", sim_horizon + 1L)
              )
          )
      )
  )

# plot the trajectory
if (show_plots) {
  get_palette <- scales::brewer_pal(type = "qual", palette = 1)
  p <- plot(patient, fill = get_palette, verbose = TRUE)
  print(p)
  # htmlwidgets::saveWidget(p, "patient_trajectory.html") # Save as html
}

# Set simulation --------------------------------------------------------------

# 1. We are only interested about people that will seek healthcare resources
# From now on cases only represent those that will use healthcare resources at some point
infected[prb_seeking_healthcare, on = "agegroup", cases := cases * i.prb]

# convert cases to integers with special treatment for <1 case
infected[cases < 1L, cases := rbinom(.N, 1L, cases)]
infected[, cases := as.integer(round(cases))] # NOTE: I assume cases < int.max

# dataset for validation with all the cases to run through the DES pathway
in_val <- infected[, .(cases_in = sum(cases)), keyby = .(loc, iteration)]


# 2. Estimate when they will seek healthcare resources. People will usually not
# do that immediately after they develop symptoms.
# First, we convert date to continuous with 0 = min date. Time is assumed in
# days. I.e. 1 is a day and 1/24 = 0.042 is an hour
starting_date <- infected[, min(date)]
infected[, time := date - ..starting_date]

# Then we need to add up the time until cases seek healthcare. For that we
# assume a Gamma distr. Here things get tricky, as we need to convert this into
# a microsimulation. The following code will create one row for every case, and
# then will add in the time until healthcare seeking. If the number of
# iterations is high, then this will produce a very long table and tables with
# more than 2 billions of rows are not allowed. Hence from now on I will loop
# over the iterations and will produce a table per iteration
infected <- split(infected, by = "loc")
infected <- lapply(infected, function(x) x[rep(1:.N, cases)][, cases := NULL])
invisible(lapply(infected, function(x) x[, time :=
                                           time + time_to_seek_healthcare(.N)]))
invisible(lapply(infected, setkeyv, "time"))

# Run simulation --------------------------------------------------------------
# Remember to delete files in the Output folder!!!
if (delete_previous_results) file.remove(list.files(normalizePath("./Outputs/"), full.names = TRUE))

# For Windows machines only use 1 core
if (Sys.info()[1] == "Windows") n_cores <- 1L

# Nested loop over locations (outer) and iterations (inner, parallelised)

for (iter_ in names(infected)) {
  print(iter_)
  hospital <- simmer("Hospital", log_level = 0L)
  out <- mclapply(infected[[iter_]][, sort(unique(iteration))], function(i) {
    hospital %>%
      add_resource("A&E", 1e7) %>% # Almost infinite capacity
      add_resource("A&E_discharge", 1e7) %>%
      add_resource("bed1", 1e7) %>% # Almost infinite capacity
      add_resource("bed1_discharge", 1e7) %>%
      add_resource("ICU", 1e7) %>% # Almost infinite capacity
      add_resource("ICU_discharge", 1e7) %>%
      add_resource("in-hospital_death", 1e7) %>%
      add_dataframe(
        "patient",
        patient,
        infected[[iter_]][iteration == i],
        time = "absolute",
        col_attributes = "agegroup",
        mon = 2L,
        batch = 50
      ) %>%
      run(sim_horizon) %>% # , progress=progress::progress_bar$new()$update
      wrap()
  },
  mc.preschedule = FALSE, # !!!If TRUE, produces wrong results
  mc.cores = n_cores)

  # Collect resources --------------------------------------------------------------
  resource_utilisation <- get_mon_resources(out)
  setDT(resource_utilisation)
  resource_utilisation[, time_int := as.IDate(floor(time) + starting_date)]
  setnames(resource_utilisation, "replication", "iteration")
  # resource_utilisation[area_lookup, on = "code", loc := i.loc]
  # resource_utilisation[, code := NULL]

  fwrite(
    resource_utilisation[, .(server = max(server), loc = iter_),
                         keyby = .(resource, time_int, iteration)],
    normalizePath("./Outputs/resource_utilisation.csv"),
    append = TRUE
  )

  # Collect individual trajectories -----
  if (collect_indiv_traj) {
    patient_trajectories <- get_mon_arrivals(out, TRUE, TRUE)
    setDT(patient_trajectories,
          key = c("start_time", "name", "replication"))
    setnames(patient_trajectories, "replication", "iteration")

    patient_trajectories <-
      unique(patient_trajectories,
             by = c("name", "start_time", "resource", "iteration"))
    patient_trajectories[, loc := iter_]
    fwrite(
      patient_trajectories,
      normalizePath("./Outputs/patient_trajectories.csv"),
      append = TRUE
    )
  }
}



# Collect resources --------------------------------------------------------------
resource_utilisation <-
  fread(
    normalizePath("./Outputs/resource_utilisation.csv"),
    stringsAsFactors = FALSE,
    colClasses = list(
      IDate = "time_int",
      character = "resource",
      numeric = "server",
      factor = "loc"
    ),
    key = c("resource", "time_int", "loc")
  )

# Expand with missing dates
tt <- CJ(time_int = as.IDate(starting_date:(starting_date + sim_horizon)),
         resource = unique(resource_utilisation$resource),
         iteration = unique(resource_utilisation$iteration),
         loc = unique(resource_utilisation$loc))
resource_utilisation <- resource_utilisation[tt, on = .NATURAL, ]

# Fill server with last obs carried forward (locf)
setkey(resource_utilisation, time_int, resource, loc, iteration)
resource_utilisation[, server := nafill(server, "locf"), by = .(resource, loc, iteration)]
# Set to 0 rows before the first recorded event
setnafill(resource_utilisation, "c", 0L, cols = "server", )


# Internal validation --------------------------------------------------------
# check that all cases that went into the model (in_val) got out. It only works
# when the simulation horizon is longer that the input data so all cases are
# accumulated in the "absorbing states". Otherwise, double counting occurs
# because the output in the non-absorbing states is discretised and do not
# reflect the individual trajectories
out_val <- resource_utilisation[time_int == max(time_int), .(cases_out = sum(server)), keyby = .(loc, iteration)]

val <- in_val[out_val, on = .(loc, iteration)]
val[, error := (cases_out - cases_in)/cases_in]
val[, summary(error)]
if (min(val$error) < 0 || max(val$error) > 0.1) stop("Relative error between cases in the input and output exceeds 10% or cases in the output is less than cases in the input. Note that some error is expected with output cases being slightly more than input cases because some double counting exists. I.e. when patients move from ICU to level 1 bed they may be counted twice.")


# remove dates after the last date of the input file
infected[[iter_]][, max(date)]
resource_utilisation <- resource_utilisation[time_int <= infected[[iter_]][, max(date)]]

# Summarise iterations in median with uncertainty intervals ------------------
resource_utilisation <-
  resource_utilisation[, as.list(round(fquantile(server, probabilities))),
                       keyby = c("time_int", "loc", "resource")]


# beautify col names
setnames(
  resource_utilisation,
  c("time_int", "loc", "resource", paste0("V", seq_len(
    length(probabilities)
  ))),
  c(
    "Date",
    "Location",
    "Resource",
    "Median",
    paste0(percent(probabilities[-1]), " UI")
  )
)

resource_utilisation[, Resource := factor(Resource)]
fwrite(resource_utilisation, normalizePath("./Outputs/summarised_output.csv"))

# Plots -------------------------------------------------------------------
ggplot(
  resource_utilisation[Resource %in% c("bed1", "ICU", "in-hospital_death")],
  aes(
    x = Date,
    y = Median,
    col = Resource,
    fill = Resource,
    ymin = `10.0% UI`,
    ymax = `90.0% UI`
  )
) +
  geom_line() +
  geom_ribbon(alpha = 1/4, linetype = 0) +
  facet_wrap( ~ Location)


ggplot(
  resource_utilisation[Resource %in% c("A&E_discharge", "ICU_discharge", "bed1_discharge")],
  aes(
    x = Date,
    y = Median,
    col = Resource,
    fill = Resource,
    ymin = `10.0% UI`,
    ymax = `90.0% UI`
  )
) +
  geom_line() +
  geom_ribbon(alpha = 1/4, linetype = 0) +
  facet_wrap( ~ Location)
