library(shiny)
library(shinyjs)
library(shinythemes)
library(shinyalert)
library(shinyBS)
library(shinyWidgets)
library(plotly)
library(DT)

# Define UI for application that draws a histogram
shinyUI(
    shiny::navbarPage("KMD Modeling Case Study",id = "page",
                      header = tagList(
                          conditionalPanel("input.page == 'setup'",
                                           fluidRow(
                                               column(10, offset = 1,
                                                      radioGroupButtons("rdo_setup_tab",NULL,
                                                                        choices = list(
                                                                            "Chemical and Physiology"="chem_physio",
                                                                            "Metabolism Parameters"="metab",
                                                                            "Urinary Clearance Parameters"="urinary",
                                                                            "Simulation Setup"="sim"
                                                                        ),
                                                                        selected = "chem_physio",
                                                                        width = "100%"
                                                      )
                                               )
                                           )
                                           ),
                          
                          fluidRow(
                              tags$br()
                          )
                          
                      ),
                      tabPanel(title = NULL,icon = icon("home"),
                               includeHTML("www/home.html")
                               ),
                      tabPanel(title = "Setup",value ="setup",
                               tabsetPanel(id = "setup_tabs",type= "hidden",
                                           tabPanelBody(value = "chem_physio",
                                                        fluidRow(
                                                            column(3,
                                                                   selectInput("sel_chem","Select Chemical",
                                                                               choices =c("Chemical A"="ChemA","Chemical B"="ChemB",
                                                                                          "Chemical C"="ChemC","Chemical D"="ChemD",
                                                                                          "Chemical E"="ChemE"),width = "100%")
                                                            ),
                                                            column(3,
                                                                   selectInput("sel_org","Select Physiology",
                                                                               choices =c("Rat"="ra","Human"="ha"),width = "100%")
                                                            )
                                                        ),
                                                        fluidRow(
                                                          tags$h3("Chemical Properties")
                                                        ),
                                                        fluidRow(
                                                            column(4,
                                                                   numericInputIcon("num_MW","Molecular Weight",
                                                                                    icon = list("g/mol"),value = 0,
                                                                                    width= "100%")
                                                                   ),
                                                            column(4,
                                                                   numericInput("num_lkow","Log10 Ocatnol Water Partition",
                                                                                value = 1,width = "100%")
                                                                   )#,
                                                            # column(4,
                                                            #        numericInputIcon("num_MWmet","Molecular Weight of the metabolite",
                                                            #                         icon = list("g/mol"),value = 1,width = "100%")
                                                            #        )
                                                        ),
                                                        fluidRow(
                                                          tags$h3("Physiology")
                                                        ),
                                                        fluidRow(
                                                            column(4,
                                                                   numericInputIcon("num_bw","Body Weight",icon = list("kg"),
                                                                                    value = 0.35,min = 0.0001, 
                                                                                    help_text ="Weight should be greater than 0" ,width = "100%")
                                                                   ),
                                                            column(4,
                                                                   numericInputIcon("num_qcc","Cardiac Output",icon= list("L/h/ kg.BW \U00BE"),
                                                                                    value = 15.2,width = "100%")
                                                                   ),
                                                            column(4,
                                                                   numericInput("num_hct","Hematocrit Factor",value = 0.42,min = 0, max = 1, width = "100%")
                                                        )
                                                        ),
                                                        fluidRow(
                                                            column(3,
                                                                   numericInput("num_vbldc","Fractional Blood Volume",
                                                                                value= 0.074, min = 0, max =1 , width = "100%")
                                                                   ),
                                                            column(3,
                                                                   numericInput("num_vlivc","Fractional Liver Volume",
                                                                                value= 0.0387, min = 0, max =1 , width = "100%")
                                                            ),
                                                            column(3,
                                                                   numericInput("num_vrpfc","Fractional Rapidly Perfused Tissue Volume",
                                                                                value= 0.0647, min = 0, max =1 , width = "100%")
                                                            ),
                                                            column(3,
                                                                   numericInput("num_vspfc","Fractional Slowly Perfused Tissue Volume",
                                                                                value= 0.6925, min = 0, max =1 , width = "100%")
                                                            )
                                                        ),
                                                        fluidRow(
                                                            column(3,
                                                                   numericInput("num_qlivc","Fractional Liver Blood Flow",
                                                                                value= 0.183, min = 0, max =1 , width = "100%")
                                                            ),
                                                            column(3,
                                                                   numericInput("num_qrpfc","Fractional Rapidly Perfused Tissue Blood Flow",
                                                                                value= 0.58, min = 0, max =1 , width = "100%")
                                                            ),
                                                            column(3,
                                                                   numericInput("num_qspfc","Fractional Slowly Perfused Tissue Blood Flow",
                                                                                value= 0.23, min = 0, max =1 , width = "100%")
                                                            ),
                                                            column(3,
                                                                   numericInputIcon("num_vurinec","Urine Production",icon = list("L/kg BW/day"),
                                                                                value= 0.012, min = 0, width = "100%")
                                                            )
                                                        ),
                                                        fluidRow(
                                                          tags$h3("Partitions")
                                                        ),
                                                        fluidRow(
                                                            column(4,
                                                                   numericInput("num_pliv","Liver Parition Coefficient",
                                                                                value= 0.183, min = 0, max =1 , width = "100%")
                                                            ),
                                                            column(4,
                                                                   numericInput("num_prpf","Rapidly Perfused Tissue Partition Coefficient",
                                                                                value= 0.58, min = 0, max =1 , width = "100%")
                                                            ),
                                                            column(4,
                                                                   numericInput("num_pspf","Slowly Perfused Tissue Parition Coefficient",
                                                                                value= 0.23, min = 0, max =1 , width = "100%")
                                                            )
                                                        )
                                                        ),
                                           tabPanelBody(value = "metab",
                                                        fluidRow(
                                                            column(4,
                                                                   numericInput("num_vmax","Maximum Metabolic Rate (\U00B5mol/h)",
                                                                                value = 0,width = "100%")
                                                            ),
                                                            column(4,
                                                                   numericInput("num_km","Michaelis Menten Constant (\U00B5M)",
                                                                                value = 1,width = "100%")
                                                                   
                                                            ),
                                                            column(4,
                                                                   numericInput("num_vkm1","Hepatic Clearence (L/h)",
                                                                                value = 0,width = "100%")
                                                            )
                                                        )
                                                        ),
                                           tabPanelBody(value = "urinary",
                                                        fluidRow(
                                                          column(4,
                                                                 numericInput("num_vmaxu","Maximum Urinary Excretion Rate (\U00B5mol/h)",
                                                                              value = 0,width = "100%")
                                                          ),
                                                          column(4,
                                                                 numericInput("num_kmu","Michaelis Menten Constant (\U00B5M)",
                                                                              value = 1,width = "100%")
                                                                 
                                                          ),
                                                          column(4,
                                                                 numericInput("num_vke1","First-order Urinary Clearance Rate (L/h)",
                                                                              value = 0,width = "100%")
                                                          )
                                                        )
                                                        
                                                        ),
                                           tabPanelBody(value = "sim",
                                                        fluidRow(
                                                          column(3,
                                                                 selectInput("sel_testExpoRoute","Simulation Exposure Route",
                                                                             choices = c("IV infusion"="iv","Oral Bolus Dose"="boral"),
                                                                             width = "100%"
                                                                 )
                                                          ),
                                                          conditionalPanel("input.sel_testExpoRoute=='iv'",
                                                                           
                                                                           column(3, 
                                                                                  numericInputIcon("num_testIVDose",
                                                                                                   "IV Exposure for Time Course Simulation",icon = list("mg/h"),
                                                                                                   value = 10,
                                                                                                   width = "100%")
                                                                           ),
                                                                           column(3,
                                                                                  numericInputIcon("num_testIVLen",
                                                                                                   "Length of IV Infusion",
                                                                                                   icon = list("h/day"),
                                                                                                   value = 24,
                                                                                                   min = 0.1, 
                                                                                                   max = 24,
                                                                                                   help_text = "Value should be between 0.1 and 24",
                                                                                                   width = "100%"
                                                                                  )
                                                                           )
                                                                           
                                                          ),
                                                          conditionalPanel("input.sel_testExpoRoute == 'boral'",
                                                                           
                                                                           column(2,
                                                                                  numericInputIcon("num_testOralDose", 
                                                                                                   "Bolus Exposure for Time Course Simulation",
                                                                                                   icon = list("mg/kg BW/day"),
                                                                                                   value = 10,
                                                                                                   width= "100%")
                                                                           ),
                                                                           column(2,
                                                                                  numericInputIcon("num_ka","Oral Absortion Rate",icon= list("/h"),
                                                                                                   value = 5,min = 0,help_text = "Value should be greater than 0"
                                                                                  )
                                                                           ),
                                                                           column(2,
                                                                                  numericInputIcon("num_fa",
                                                                                                   "Fraction Abrorbed", 
                                                                                                   value = 1, min = 0, max = 1)
                                                                                  )
                                                                           
                                                          )
                                                        ),
                                                        fluidRow(
                                                          column(3,
                                                                 selectInput("sel_metabType","Select Metabolism Type",
                                                                             c("Saturable"="sat","Linear"="lin"),
                                                                             width = "100%"
                                                                 )
                                                          ),
                                                          column(3,
                                                                 selectInput("sel_uclType","Select Urinary Clearance Type",
                                                                             c("Saturable"="sat",
                                                                               "Linear First Order"="lin"),width = "100%"
                                                                 )
                                                          ),
                                                          column(3,
                                                                 numericInputIcon("num_testSimDuration","Simulation Duration",
                                                                                  icon = list("h"),
                                                                                  value= 2160, 
                                                                                  min = 0.1,
                                                                                  help_text = "Value has to be greater than 0",
                                                                                  width = "100%"
                                                                 )
                                                          )
                                                        ),
                                                        fluidRow(
                                                          column(6,
                                                                 numericRangeInput("numrange_expo",
                                                                                   "Select Exposure Range for Dose Responses Plots",
                                                                                   value = c(0.1,10),separator = "to",
                                                                                   width = "100%"
                                                                                   
                                                                                   
                                                                 )
                                                          ),
                                                          column(6,
                                                                 numericInput("num_numexpos","Number of Exposures",
                                                                              value = 50,width = "100%")
                                                          )
                                                        ),
                                                        
                                                        fluidRow(
                                                          column(4,offset=2,
                                                                 actionButton("btn_runTestSim","Run Time Course Simulation",width = "100%")
                                                          ),
                                                          column(4, offet = 2, 
                                                                 actionButton("btn_runDRSim","Run Dose Reponse Simulation",width = "100%")
                                                          )
                                                        ),
                                                        bsCollapse(id="sim_panels",
                                                                   bsCollapsePanel("Time Course Plot",value = "test_sims",
                                                                                   fluidRow(
                                                                                     column(10,offset = 1,
                                                                                            plotlyOutput("plt_testSim")
                                                                                            )
                                                                                       )
                                                                                   ),
                                                                   bsCollapsePanel("Dose Response Simulation Plots", value = "DR_sims",
                                                                                   tabsetPanel(
                                                                                     tabPanel("Plot",
                                                                                              fluidRow(
                                                                                                column(1,
                                                                                                       selectInput("sel_yaxis",NULL,
                                                                                                                   choices = c("log","linear"))
                                                                                                ),
                                                                                                column(10,offset = 1,
                                                                                                       plotlyOutput("plt_DRSim")
                                                                                                )
                                                                                                
                                                                                              )
                                                                                              ),
                                                                                     tabPanel("Data",
                                                                                              fluidRow(
                                                                                                column(10,offset = 1,
                                                                                                       DTOutput("tble_DRSim")
                                                                                                       )
                                                                                                )
                                                                                              )
                                                                                   )
                                                                                   
                                                                                   )
                                                                       
                                                         
                                                            
                                                        )
                                                        )
                                           )
                               ),
                      # tabPanel("Output", value = "output"),
                      tabPanel(NULL,value = "off",icon = icon("power-off"))
                      )
    )

