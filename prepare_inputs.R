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


# Read output of Lancaster transmission model --------------------------------------
infected <- fread(normalizePath("./Inputs/infections.csv"),
            stringsAsFactors = FALSE,
            colClasses = list(
              IDate = "date",
              character = "loc",
              numeric = "cases"
            ))


# Plot input file
if (show_plots) {
setkey(infected, agegroup)
tt <- infected[, as.list(round(fquantile(cases, probabilities))),
                     by = c("date", "loc", "agegroup")]
# beautify col names
setnames(
  tt,
  c("date", "loc", "agegroup", paste0("V", seq_len(
    length(probabilities)
  ))),
  c(
    "Date",
    "Location",
    "Agegroup",
    "Median",
    paste0(percent(probabilities[-1]), " UI")
  )
)

ggplot(
  tt,
  aes(
    x = Date,
    y = Median,
    col = Agegroup,
    fill = Agegroup,
    ymin = `20.0% UI`,
    ymax = `80.0% UI`
  )
) +
  geom_line() +
  geom_ribbon(alpha = 1/4, linetype = 0) +
  facet_grid(Agegroup ~ Location, scales = "fixed")
}


# I need to replace loc with an int, to be used later in parallelisation. I
# create a look-up table to convert back to area codes for outputs
# area_lookup <- infected[, .(loc = sort(unique(loc)))][, code := seq_len(.N)]
# infected[area_lookup, on = "loc", loc := i.code]


# Prepare user inputs --------------------------------------
prb_seeking_healthcare <- fread(normalizePath("./Inputs/prb_seeking_healthcare.csv")
)[, agegroup := as.integer(gsub("*\\-.*|*\\+.*", "", agegroup))]
stopifnot(all(prb_seeking_healthcare$prb >= 0)) # error checking
if(any(prb_seeking_healthcare$prb > 1)) message("prb > 1 in healthcare seeking")


prb_hospitalisation <- fread(normalizePath("./Inputs/prb_hospitalisation.csv")
)[, agegroup := as.integer(sub("*\\-.*|*\\+.*", "", agegroup))]
# Change denominator from all cases to those seeking healthcare resources
prb_hospitalisation[prb_seeking_healthcare, on = "agegroup", prb := prb/i.prb]
stopifnot(all(prb_hospitalisation$prb >= 0)) # error checking
if(any(prb_hospitalisation$prb > 1)) message("prb > 1 in hospitalisation")

prb_hospitalisation_vec <- prb_hospitalisation$prb
names(prb_hospitalisation_vec) <- prb_hospitalisation$agegroup

get_traj_hospitalisation <- function(agegrp, dt = prb_hospitalisation_vec) {
  agegrp <- as.character(agegrp)
  if (length(agegrp) != 1L) stop()
  if (!agegrp %in% names(dt)) stop()
  rbinom(1L, 1L, dt[[agegrp]]) + 1L # 1 = discharge, 2 = hospitalise
}


prb_icu <- fread(normalizePath("./Inputs/prb_icu.csv")
)[, agegroup := as.integer(sub("*\\-.*|*\\+.*", "", agegroup))]
# Change denominator from all cases to those hospitalised
tt <- fread(normalizePath("./Inputs/prb_hospitalisation.csv")
)[, agegroup := as.integer(sub("*\\-.*|*\\+.*", "", agegroup))]
prb_icu[tt, on = "agegroup", prb := prb/i.prb]
stopifnot(all(prb_icu$prb >= 0)) # error checking
if(any(prb_icu$prb > 1)) message("prb > 1 in ICU")
rm(tt)

prb_icu_vec <- prb_icu$prb
names(prb_icu_vec) <- prb_icu$agegroup

get_traj_icu <- function(agegrp, dt = prb_icu_vec) {
  agegrp <- as.character(agegrp)
  if (length(agegrp) != 1L) stop()
  if (!agegrp %in% names(dt)) stop()
  rbinom(1L, 1L, dt[[agegrp]]) + 1L # 1 = normal bed, 2 = ICU
}


prb_icu_death <- fread(normalizePath("./Inputs/prb_icu_death.csv")
)[, agegroup := as.integer(sub("*\\-.*|*\\+.*", "", agegroup))]
stopifnot(all(prb_icu_death$prb >= 0)) # error checking
if(any(prb_icu_death$prb > 1)) message("prb > 1 in hospitalisation")

prb_icu_death_vec <- prb_icu_death$prb
names(prb_icu_death_vec) <- prb_icu_death$agegroup

get_traj_icu_death <- function(agegrp, dt = prb_icu_death_vec) {
  agegrp <- as.character(agegrp)
  if (length(agegrp) != 1L) stop()
  if (!agegrp %in% names(dt)) stop()
  rbinom(1L, 1L, 1 - dt[[agegrp]]) + 1L # 1 = in-hospital death, 2 = de-escalate
}


prb_nonicu_death <- fread(normalizePath("./Inputs/prb_nonicu_death.csv")
)[, agegroup := as.integer(sub("*\\-.*|*\\+.*", "", agegroup))]
stopifnot(all(prb_nonicu_death$prb >= 0)) # error checking
if(any(prb_nonicu_death$prb > 1)) message("prb > 1 in hospitalisation")

prb_nonicu_death_vec <- prb_nonicu_death$prb
names(prb_nonicu_death_vec) <- prb_nonicu_death$agegroup

# get_traj_nonicu_death <- function(agegrp, dt = prb_nonicu_death_vec) {
#   agegrp <- as.character(agegrp)
#   if (length(agegrp) != 1L) stop()
#   if (!agegrp %in% names(dt)) stop()
#   rbinom(1L, 1L, dt[[agegrp]]) + 1L # 1 = discharge, 2 = in-hospital death
# }

get_traj_nonicu_death <- function(agegrp, dt = prb_nonicu_death_vec) {
1L # Everyone discharged to community post ICU. Assume no in-hospital death in
   # the short period of deescalation to level 1 bed. If this assumption is 3
   # unacceptable, use the commented out function above
}
