output$reproducibilitySteps <- renderUI({
  req(reactiveVals$se)
  if(!is.null(metadata(reactiveVals$se)$steps)){
    steps <- tags$ol(
        lapply(metadata(reactiveVals$se)$steps, function(x) {
          tags$li(x)
        })
      )
    div(steps,  style = "font-size: large;")
    
  } else {
      div("No Steps Executed With the Current Data !", style = "font-size: large;")
  }
})