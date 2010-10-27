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
max.n.observed <- 9
observed.names = c("parent", paste("M", 1:(max.n.observed - 1), sep=""))

prg <- ggroup(horizontal=FALSE, cont = pr)
prl <- glayout(cont = prg)
prl[1,1] <- glabel("Number of observed variables", cont=prl)
prl[1,2] <- (n.observed.gw = gcombobox(
  1:max.n.observed, 
  handler = function(h, ...) {
    n.observed <- svalue(n.observed.gw)
    filter.expression <- paste("^[1-", n.observed, "]$", sep="")
    observed.gw$filter("Index", filter.expression)
  },
  cont=prl))

observed.gw <- gdf(
  items = data.frame(Index = 1:max.n.observed, 
    Name = observed.names, stringsAsFactors=FALSE),
  name = "Names of observed variables",
  width = 500, height=250, cont=prg)

addSpring(prg)
addSpace(prg, 15)
pr.ok <- gbutton("ok", 
  handler=function(h,...) {
      galert(getwd())
    }, container=prg)

visible(pr) <- TRUE

# Dataset definition
ds <- gexpandgroup("Datasets", cont=w)
s <- gdf(data.frame( Author = c("Meier"), Year = "3", Title = "4"),
  name = "Studies",
cont=ds)

size(s) <- c(800,200)
visible(ds) <- TRUE

# Model definition
ms <- gexpandgroup("Models", cont=w)
visible(ms) <- FALSE

# Evaluation window
mw <- gframe("Evaluations", cont=w)

# Status bar
gstatusbar("Powered by gWidgetsWWW", cont = w)
visible(w) <- TRUE 