library(caret)
library(DT)
library(dplyr)
library(ggplot2)
library(plotly)
library(randomForest)
library(rsconnect)
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(tidyverse)

data <- read.csv("49genes.csv")
df <- data %>% gather(key = "ID", value = "value", -X)
model <- readRDS("rf_fit_final_model.rds")

ui <- navbarPage("Kidney Pro",
                 theme = shinytheme("cerulean"),
                 tabPanel("Home",
                          h3("Welcome to Kidney Pro!"),
                          p("This app provides an overview of your patients' information and insights into the
                             rejection status of your patient along with treatment suggestions based on the rejection type.
                             Please upload a CSV file containing your patients' gene expression data.
                             This file should be in the format of gene symbols as rows and samples as columns."),
                          h3("Target Audience"),
                          p("A junior clinician nephrologist who is going to test biopsy samples and provide 
                            a detailed treatment and care plan for a patient who has recently undergone a kidney transplant operation.")
                 ),
                 tabPanel("Prediction",
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("file", "Select a CSV File",
                                        multiple = FALSE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              selectInput("id", label = "Select Patient ID", choices = NULL),
                              actionButton("action", label = "Predict")
                            ),
                            mainPanel(
                              uiOutput("text")
                            )
                          )         
                 ),
                 tabPanel("Insights",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("id2", label = "Select Patient ID", choices = NULL),
                              pickerInput("gene", "Select Gene(s):",
                                          choices = sort(unique(df$X)), selected = c("PLA1A", "BATF", "GBP5", "CCL4", "ROBO4", "CXCR6"), options = list(`actions-box` = TRUE), multiple = TRUE)
                            ),
                            mainPanel(
                              plotlyOutput("plot"),
                              p("Nearer to mean value: Higher likelihood of a Stable diagnosis. Further from mean value (outlier): Indicative of some type of Rejection.")
                            )
                          )
                 )
)

server <- function(input, output, session) {
  observeEvent(input$file, {
    req(input$file)
    test <- read.csv(input$file$datapath, row.names = 1)
    updateSelectInput(session, "id", choices = colnames(test), selected = "GSM3701773_2")
    updateSelectInput(session, "id2", choices = colnames(test), selected = "GSM3701773_2")
  })
  
  output$plot <- renderPlotly({
    req(input$id2)
    req(input$gene)
    
    if (length(input$gene) == 0) {
      return(NULL)
    }
    
    df2 <- df %>% filter(ID %in% input$id2, X %in% input$gene)
    df %>%
      filter(X %in% input$gene) %>%
      ggplot() + geom_boxplot(aes(x = X, y = value)) +
      geom_point(data = df2, aes(x = X, y = value), color = "red", size = 2) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Gene", y = "",
           title = "Gene Expression of Samples vs Normal Gene Distribution")
  })
  
  output$text <- renderUI({
    input$action
    req(input$file)
    test <- read.csv(input$file$datapath, row.names = 1)
    features <- rownames(varImp(model))
    X_test_final_model <- test[rownames(test) %in% features, ]
    X_test_final_model <- t(X_test_final_model)
    preds <- predict(model, X_test_final_model)
    probs <- predict(model, X_test_final_model, type = "prob")
    df <- data.frame(gsm = colnames(test), predicted = preds)
    df <- cbind(df, probs)
    pred <- df$predicted[df$gsm == input$id]
    if (length(pred) == 0) {
      tags$div( 
        h4("Not found.")
      )
    } else {
      if (pred == "ABMR") {
        prob <- df$ABMR[df$gsm == input$id]
        tags$div(
          h4(paste0("The predicted class is: ", pred, " (", round(prob * 100, 2), "%)"), style = "color:red"),
          hr(),
          h4("ABMR"),
          h5("Definition"),
          tags$ul(
            tags$li("ABMR, also known as antibody-mediated rejection, it is a B cell mediated process characterised by
                     the production of IgG antibodies against the transplanted kidney. This type of rejection is one
                     of the main barriers to the long-term survival of kidney transplant patients."),
            tags$li("The morphologic nature of endothelial cell injury in acute ABMR demonstrates platelet
                     aggregation, thrombotic microangiopathy (TMA) and neutrophilic accumulation, resulting in
                     an early pattern of cellular necrosis and a relatively rapid decline in allograft function.")
          ),
          h5("Types"),
          tags$ul(
            tags$li("Acute/active ABMR: Characterized by rapid onset and involves microvascular inflammation,
                     evidence of current or recent antibody interaction with vascular endothelium, and serologic
                     evidence of DSAs."),
            tags$li("Chronic active ABMR: This results from a repetitive pattern of chronic thrombotic events and
                     inflammatory changes, which result in cellular injury and repair and can be diagnosed from
                     C4d (a degradation product of the complement pathway that binds covalently to the
                     endothelium) and DSA (anti-HLA donor-specific antibodies).")
          ),
          h5("Diagnosis"),
          tags$ul(
            tags$li("Histologic evidence of tissue injury."),
            tags$li("Immunohistochemical staining (e.g., C4d staining in peritubular capillaries) and other
                     markers of endothelial damage."),
            tags$li("Molecular diagnostics such as increased expression of endothelial activation and injury
                     transcripts (ENDATs).")
          ),
          h5("Treatment"),
          tags$ul(
            tags$li("Immunosuppressive therapies to control the immune response."),
            tags$li("Plasmapheresis to remove circulating antibodies."),
            tags$li("Intravenous immunoglobulin (IVIG) to neutralize circulating antibodies."),
            tags$li("Anti-B cell therapies (like rituximab) to deplete B cells that produce antibodies."),
            tags$li("Anti-thymocyte globulin (ATG) for more severe cases."),
            tags$li("Corticosteroids to reduce inflammation and immune response.")
          ),
          em("Note: these predictions and recommendations are for reference only.", style = 'color:red')
        )
      } else if (pred == "TCMR") {
        prob <- df$TCMR[df$gsm == input$id]
        tags$div(
          h4(paste0("The predicted class is: ", pred, " (", round(prob * 100, 2), "%)"), style = "color:red"),
          hr(),
          h4("TCMR"),
          h5("Definition"),
          tags$ul(
            tags$li(
              "TCMR, also known as T cell-mediated rejection, refers to the rejection due to a T cell-mediated immune response. It is characterised
               by the infiltration of T cells into the graft, causing damage and dysfunction."
            )
          ),
          h5("Types"),
          tags$ul(
            tags$li("There are no explicitly categorised types of TCMR, but different severities of rejection are classified with the
                     Banff classification, such as Suspicious (Borderline) for Acute TCMR, Acute TCMR IA, IB, IIA, and
                     higher. These categories indicate varying degrees of T cell infiltration and tissue damage.")
          ),
          h5("Diagnosis"),
          tags$ul(
            tags$li("Histological assessment using biopsy to identify T cell infiltration."),
            tags$li("Grading of rejection severity using the Banff classification system, which can identify
                     different grades of rejection based on tissue samples."),
            tags$li("Monitoring graft function through clinical markers such as serum creatinine levels, though
                     these are considered less sensitive.")
          ),
          h5("Treatment"),
          tags$ul(
            tags$li("High-dose corticosteroids to suppress the immune response."),
            tags$li("Increased maintenance immunosuppression with medications such as tacrolimus and
                     mycophenolic acid."),
            tags$li("Antibody therapies like thyroglobulin for severe cases, targeting T cells directly."),
            tags$li("Adjusting immunosuppressive regimens to prevent further episodes of rejection.")
          ),
          em("Note: these predictions and recommendations are for reference only.", style = 'color:red')
        )
      } else {
        prob <- df$Stable[df$gsm == input$id]
        tags$div(
          h4(paste0("The predicted class is: ", pred, " (", round(prob * 100, 2), "%)"), style = "color:red"),
          hr(),
          em("Note: this prediction is for reference only.", style = 'color:red')
        )
      }
    }
  })
}

# Create Shiny app
shinyApp(ui = ui, server = server)