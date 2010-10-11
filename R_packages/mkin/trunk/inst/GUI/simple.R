# Work in progress - simple GUI for mkin
w <- gwindow('Simple R GUI for kinetic evaluations', visible=FALSE) 

# Handlers for the Project menu
new.project.handler = function(h,...) gmessage("Called new project handler", parent=w)
open.project.handler = function(h,...) gmessage("Called open project handler", parent=w)

# Actions
actionlist = list(
  new = gaction(label="New", icon="new", 
    handler = new.project.handler, parent = w),
  open = gaction(label="Open", icon="open", 
    handler = open.project.handler, parent = w)
)

# The menu bar
menulist <- list(Project = list(
  new = actionlist$new,
  open = actionlist$open
))

mb = gmenu(menulist, cont=w)

# Project definition
pr <- gexpandgroup("Project", cont=w)

n.observed <- 1
max.n.observed <- 20
observed.names = c("parent", paste("M", 1:19, sep=""))

prg <- ggroup(horizontal=FALSE, cont = pr)
prl <- glayout(cont = prg)
prl[1,1] <- glabel("Number of observed variables", cont=prl)
prl[1,2] <- (n.observed.gw = gcombobox(
  1:max.n.observed, 
  handler = function(h, ...) {
    n.observed <- svalue(n.observed.gw)
  },
  cont=prl))

observed.gw <- gdf(
  items = data.frame(Index = 1:max.n.observed, Name = observed.names, stringsAsFactors=FALSE),
  name = "Names of observed variables",
  width = 500, height=500, cont=prg)
visible(observed.gw) <- FALSE

visible(pr) <- TRUE

# Dataset definition
ds <- gexpandgroup("Datasets", cont=w)
visible(ds) <- FALSE

# Model definition
ms <- gexpandgroup("Models", cont=w)
visible(ms) <- FALSE

# Evaluation window
mw <- gframe("Evaluations", cont=w)

visible(w) <- TRUE 
