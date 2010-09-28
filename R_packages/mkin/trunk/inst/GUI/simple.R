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

# Layout
ds <- gexpandgroup("Datasets", cont=w)
visible(ds) <- FALSE
ms <- gexpandgroup("Models", cont=w)
visible(ds) <- FALSE
mw <- gframe("Main window", cont=w)

visible(w) <- TRUE 
