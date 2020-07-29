#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinybusy)
library(NonCompart)
library(RSQLite)
library(deSolve)
library(pracma)
library(magrittr)
library(plotly)
source("R/KMD_model_event_inits.R")
source("R/QSAR Models.R")
source("R/Parameters.R")

# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
    #create database connection
    conn <- dbConnect(SQLite(),"Database/ModelDb.sqlite")
    #Identify the model to run
    C_mName <- file.path("Model","KMD_model_event.c")
    system(paste("R CMD SHLIB ", C_mName, sep = ""))
    dll_mName <- file.path("Model",
                           paste("KMD_model_event",.Platform$dynlib.ext,sep="")
    )
    
    results <- reactiveValues(pbpk=NULL,dr=NULL)
    observeEvent(input$sel_chem,{
        chem <- input$sel_chem
        vliv <- input$num_vlivc*input$num_bw
        chem_params <- getChemicalParams(chem)
        updateNumericInputIcon(session,"num_MW",value= chem_params[["MW"]])
        updateNumericInput(session,"num_lkow",value = chem_params[["lkow"]])
        updateNumericInputIcon(session,"num_vmax",value = chem_params[["vmaxc"]]*vliv)
        updateNumericInputIcon(session,"num_km",value = chem_params[["km"]])
        updateNumericInputIcon(session,"num_vkm1",value = chem_params[["vmaxc"]]*vliv/chem_params[["km"]])
        
        updateNumericInputIcon(session,"num_vmaxu",value = chem_params[["vmaxc"]]*vliv)
        updateNumericInputIcon(session,"num_kmu",value = chem_params[["km"]])
        updateNumericInputIcon(session,"num_vke1",value = chem_params[["vmaxc"]]*vliv/chem_params[["km"]])
        #update Partitions
        org <- input$sel_org
        part_coeffs <- calculatePartitions(chem_params,org)
        lapply(names(part_coeffs),function(x,coeffs){
            updateNumericInput(session,paste0("num_",x),value = coeffs[[x]])
        },part_coeffs)
    })
    observeEvent(input$sel_org,{
        org <- input$sel_org
        query <- sprintf("Select var,val from Physiology where org = '%s'",org)
        res <- dbSendQuery(conn,query)
        resDF <- dbFetch(res,n=-1)
        dbClearResult(res)
        resDF <- setNames(resDF$val,resDF$var)
        lapply(names(resDF),function(x,vals){
            updateNumericInputIcon(session,paste0("num_",x),value = vals[[x]])
        },resDF)
        #Update Paritions
        chem <- input$sel_chem
        chem_params <- getChemicalParams(chem)
        part_coeffs <- calculatePartitions(chem_params,org)
        lapply(names(part_coeffs),function(x,coeffs){
            updateNumericInput(session,paste0("num_",x),value = coeffs[[x]])
        },part_coeffs)
    })
    observeEvent({input$num_bw
        input$num_vlivc},{
            chem <- input$sel_chem
            vliv <- input$num_vlivc*input$num_bw
            chem_params <- getChemicalParams(chem)
            updateNumericInputIcon(session,"num_vmax",value = chem_params[["vmaxc"]]*vliv)
            updateNumericInputIcon(session,"num_km",value = chem_params[["km"]])
            updateNumericInputIcon(session,"num_vkm1",value = chem_params[["vmaxc"]]*vliv/chem_params[["km"]])
            
            updateNumericInputIcon(session,"num_vmaxu",value = chem_params[["vmaxc"]]*vliv)
            updateNumericInputIcon(session,"num_kmu",value = chem_params[["km"]])
            updateNumericInputIcon(session,"num_vke1",value = chem_params[["vmaxc"]]*vliv/chem_params[["km"]])
        })
    observeEvent(input$btn_runTestSim,{
        tstart <- 0
        tstop <- input$num_testSimDuration
        times <- seq(0,tstop,0.1)
        hep_metab_type <- input$sel_metabType
        ucl_type <- input$sel_uclType
        expo_route <- input$sel_testExpoRoute
        inputs <- reactiveValuesToList(input)
        names(inputs)<- sub("num_","",names(inputs))
        
        scaled_initial_values <- getModelParams(inputs,conn)
        # Set the parameters based on selection in UI
        #Hepatic Metabolism
        if(hep_metab_type == "sat"){
            scaled_initial_values[["vkm1"]]<-0
        }else{
            scaled_initial_values[["vmax"]]<-0
            scaled_initial_values[["km"]]<- 1
        }
        #Urinary Clearance
        if(ucl_type == "sat"){
            scaled_initial_values[["vke1"]]<- 0
        }else{
            scaled_initial_values[["vmaxu"]]<- 0
            scaled_initial_values[["kmu"]]<- 1
        }
        
        if(expo_route=="iv"){
            scaled_initial_values[["boral"]]<-0
            scaled_initial_values[["ivdose"]] <- input$num_testIVDose
            ivlen <- input$num_testIVLen
            if(ivlen == 24){
                event_times <- c(0)
            }else{
                #Number of replications of the event
                Nrep <- ceiling(max(times) / 24)
                #Find start and end times
                event_times <- rep(c(0, ivlen), Nrep) + rep(24 * (0:(Nrep - 1)), rep(2, Nrep))
            }
            
            
        }else{
            scaled_initial_values[["boral"]] <- input$num_testOralDose
            scaled_initial_values[["ivdose"]]<-0
            event_times <- head(seq(0,tstop,24),-1)
        }
        dyn.load(dll_mName)
        parms <-initParms(scaled_initial_values)
        y <- initStates(params)
        out <- ode(y,times,func= 'derivs',parms = parms,
                   dllname = "KMD_model_event",
                   initfunc = "initmod",
                   events = list(func = "event",time = event_times),
                   nout = length(Outputs),
                   outnames = Outputs)
        results$pbpk <- as.data.frame(out)
        dyn.unload(dll_mName)
    })
    
    observeEvent(input$btn_runDRSim,{
        shinybusy::show_modal_progress_line(value = 0,text = "Starting Simulation")
        tstart <- 0
        tstop <- input$num_testSimDuration
        times <- seq(0,tstop,0.1)
        hep_metab_type <- input$sel_metabType
        ucl_type <- input$sel_uclType
        expo_route <- input$sel_testExpoRoute
        inputs <- reactiveValuesToList(input)
        names(inputs)<- sub("num_","",names(inputs))
        
        scaled_initial_values <- getModelParams(inputs,conn)
        # Set the parameters based on selection in UI
        #Hepatic Metabolism
        if(hep_metab_type == "sat"){
            scaled_initial_values[["vkm1"]]<-0
        }else{
            scaled_initial_values[["vmax"]]<-0
            scaled_initial_values[["km"]]<- 1
        }
        #Urinary Clearance
        if(ucl_type == "sat"){
            scaled_initial_values[["vke1"]]<- 0
        }else{
            scaled_initial_values[["vmaxu"]]<- 0
            scaled_initial_values[["kmu"]]<- 1
        }
        #Get a vector of exposures to run dose response sim
        expo_range <- input$numrange_expo
        num_expos <- input$num_numexpos
        expo_vector <- pracma::logseq(expo_range[1],
                                      expo_range[2],
                                      num_expos)
        #setup empty vectors for dose response plots
        ramets = c(rep(NA,num_expos))
        raumets = c(rep(NA,num_expos))
        auc_prnts = c(rep(NA,num_expos))
        auc_mets = c(rep(NA,num_expos))
        cmax_prnts = c(rep(NA,num_expos))
        cmax_mets = c(rep(NA,num_expos))
        auc_tots = c(rep(NA,num_expos))
        for(idx in seq_along(expo_vector)){
            update_modal_progress(idx/num_expos,sprintf("Running Exposure %i of %i",idx,num_expos))
            each_dose <- expo_vector[[idx]]
            #Setup exposure
            if(expo_route=="iv"){
                scaled_initial_values[["boral"]]<-0
                scaled_initial_values[["ivdose"]] <- each_dose
                ivlen <- input$num_testIVLen
                if(ivlen == 24){
                    event_times <- c(0)
                }else{
                    #Number of replications of the event
                    Nrep <- ceiling(max(times) / 24)
                    #Find start and end times
                    event_times <- rep(c(0, ivlen), Nrep) + rep(24 * (0:(Nrep - 1)), rep(2, Nrep))
                }
                
                
            }else{
                scaled_initial_values[["boral"]] <- each_dose
                scaled_initial_values[["ivdose"]]<-0
                event_times <- head(seq(0,tstop,24),-1)
            }
            #run each simulation
            dyn.load(dll_mName)
            parms <-initParms(scaled_initial_values)
            y <- initStates(params)
            out <- ode(y,times,func= 'derivs',parms = parms,
                       dllname = "KMD_model_event",
                       initfunc = "initmod",
                       events = list(func = "event",time = event_times),
                       nout = length(Outputs),
                       outnames = Outputs)
            out <- as.data.frame(out)
            dyn.unload(dll_mName)
            #Last time index
            tlast24h_idx <- which(out$time==tstop-24)
            
            tlast_idx<- length(out$time)
            tlast24h_idx <- tlast_idx-(24/0.1)
            last24_time_array <- out$time[tlast24h_idx:tlast_idx]
            last24_cpls_array <- out$cpls[tlast24h_idx:tlast_idx]
            last24_cmet_array <- out$cmet[tlast24h_idx:tlast_idx]
  
            # get the values as needed by for dose response plots
            ramets[idx]<- max(out$ramet)
            raumets[idx]<- max(out$raumet)
            cmax_mets[idx]<- max(out$cmet)
            cmax_prnts[idx]<- max(out$cpls)
            auc_mets[idx]<- max(out$auc_cmet)-out$auc_cmet[tlast24h_idx]
            auc_prnts[idx]<- max(out$auc_cprnt)-out$auc_cprnt[tlast24h_idx]
            auc_tots[idx]<- max(out$auc_ctot)-out$auc_ctot[tlast24h_idx]
        }
        
        remove_modal_progress()
        results$dr <- data.frame("expos" = expo_vector,
                                 "ramet" = ramets,
                                 "raumet"=raumets,
                                 "cpls" = cmax_prnts,
                                 "cmet" = cmax_mets,
                                 "auc_prnt" = auc_prnts,
                                 "auc_met"= auc_mets,
                                 "auc_tot"=auc_tots)
    })
    
    
    
    
    tc_plt_data <- reactive({
        validate(need(results$pbpk,"Time Course Simulation not Run"))
        return(results$pbpk)
        
    })
    output$plt_testSim<- renderPlotly({
        plot_ly(data = tc_plt_data(),x = ~time)%>%
            add_trace(y= ~cpls,mode= "lines",type = "scatter",name = "Plasma Parent Concentration")%>%
            add_trace(y= ~cmet, mode = "lines",type = "scatter",name = "Plasma Metabolite Concentration")%>%
            layout(yaxis = list(title = "Concentration (\U00B5M)"),
                   xaxis = list(title = "Time (h)"),
                   legend=list(orientation = "h",y = 100)
                   )
        })
    
    dr_plt_data <- reactive({
        validate(need(results$dr,"Dose Response Simulation Not Run"))
        results$dr
        
    })
    output$plt_DRSim<- renderPlotly({
        plot_ly(data = dr_plt_data(),x = ~expos)%>%
            add_trace(y = ~ramet,type = "scatter",mode = "lines",name = "Rate of Metabolism",yaxis = "y2")%>%
            add_trace(y = ~raumet,type = "scatter",mode = "lines",name = "Rate Urinary Excretion",yaxis = "y2")%>%
            add_trace(y = ~auc_prnt, type = "scatter",mode = "markers",name = "Parent Plasma AUC")%>%
            add_trace(y= ~auc_met,type = "scatter",mode = "markers",name="Metabolite Plasma AUC")%>%
            add_trace(y= ~auc_tot,type = "scatter",mode = "markers",name="Total Chemical AUC")%>%
            layout(yaxis = list(side = "left",type= input$sel_yaxis,
                                title = "24h Concentration AUC at Steady State(\U00B5M.h)"),
                   yaxis2 = list(side = "right",automargin=TRUE,
                                 overlaying = "y",type = input$sel_yaxis,
                                 title = "Maximum Metabolism Rate (\U00B5mol/h)"),
                   xaxis = list(title= ifelse(input$sel_testExpoRoute=="iv",
                                              "IV Infusion (mg/h)",
                                              "Oral Bolus (mg/kg Bw/day)"),
                                type = "log"),
                   legend = list(orientation = 'h',
                                 y = 100)
                   )
    })
    
    output$tble_DRSim<-renderDT(DT::datatable(dr_plt_data(),rownames = F,
                                              colnames = c("Exposure","Maxmimum Metabolic Rate",
                                                           "Maximum Urinary Clearance Rate",
                                                           "Maxmimum Parent Concentration",
                                                           "Maximim Metabolite Concentration",
                                                           "24h AUC for Parent",
                                                           "24h AUC for Metabolite",
                                                           "24h AUC for Parent and Metabolite"),
                                              extensions = "Buttons",
                                              options = list(
                                                  dom = 'Bfrtip',
                                                  buttons = c('csv', 'excel')
                                              )
                                              ),server = F
                                )
    
    
    
    observeEvent(input$rdo_setup_tab,{
        
        updateTabsetPanel(session,"setup_tabs",selected = input$rdo_setup_tab)
        
    })
    observeEvent(input$page,{
        if (input$page == "off"){
            stopApp()
        }
    })

})
