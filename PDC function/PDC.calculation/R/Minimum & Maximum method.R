#' A PDC calculation Function
#'
#' This function allows you to calculate PDC in multiple measurements.
#' @param ID Patients' ID
#' @keywords PDC
#' @export
#' @examples
#' pdc_min_max_calc()


##############################################################################
# Author: Xiyue Li
# Analysis: Minimum & Maximum method
##############################################################################

pdc_min_max_calc <- function(ID = 'MRN', Data.Date = Sys.Date(),

                     # Data file: Medication Order
                     data_order = data_order,
                     Order.Date = 'Order.Date',
                     Start.Date = 'Start.Date',
                     End.Date = 'End.Date',

                     # Data file: Medication Dispensed
                     data_dispensed = data_dispensed,
                     Drug = 'Drug.Name',
                     Fill.Date = 'Fill.Date',
                     Days.Supplied = 'Days.Supplied',

                     # Data file: Patient Office Visit
                     data_visit = data_visit,
                     Visit.Date = 'Visit.Date',
                     Death.Date = 'Death.Date',
                     Followed.Days = 180,

                     # Data file: Patient Hospitalization
                     data_hosp = data_hosp,
                     Hosp.Start = 'Hosp.Start',
                     Hosp.End = 'Hosp.End',

                     Inclusion.Start = 'Inclusion.Start',
                     Inclusion.End = 'Inclusion.End',

                     # Using `doParallel` package to improve the speed of running the function
                     parallel = T, ncores = 4){

  ##############################################################################
  # Packages
  ##############################################################################

  library(dplyr)

  # Load the following packages if you need to execute loops in parallel
  if(parallel == T){
    library(parallel)
    library(foreach)
    library(doParallel)
    library(doSNOW)
    library(lubridate)
  }

  ##############################################################################
  # Data Manipulation
  ##############################################################################


  # Set the invalid number of supply days (i.e. Negative or Missing) to 0
  data_dispensed[,Days.Supplied][data_dispensed[,Days.Supplied]<0] <- 0
  data_dispensed[,Days.Supplied][is.na(data_dispensed[,Days.Supplied])] <- 0

  # Create a list of unique patient IDs regarding the visit table
  pids <- as.vector(as.matrix(unique(data_visit[,ID])))

  ##############################################################################
  # Define an inner function to calculate PDC for each patient
  ##############################################################################

  pdc_min_max_pid <- function(ID = 'MRN', patientID = 1,

                           # Data file: Medication Order
                           data_order = data_order,
                           Order.Date = 'Order.Date',
                           Start.Date = 'Start.Date',
                           End.Date = 'End.Date',

                           # Data file: Medication Dispensed
                           data_dispensed = data_dispensed,
                           Drug = 'Drug.Name',
                           Fill.Date = 'Fill.Date',
                           Days.Supplied = 'Days.Supplied',

                           # Data file: Patient Office Visit
                           data_visit = data_visit,
                           Visit.Date = 'Visit.Date',
                           Death.Date = 'Death.Date',
                           Followed.Days = 180,

                           # Data file: Patient Hospitalization
                           data_hosp = data_hosp,
                           Hosp.Start = 'Hosp.Start',
                           Hosp.End = 'Hosp.End',

                           Inclusion.Start = 'Inclusion.Start',
                           Inclusion.End = 'Inclusion.End'){

    ############################################################################
    # Extract data for the selected patient
    ############################################################################

    # Create a new object called `pid` for the selected patient ID
    pid <- pids[patientID]

    # Medication Order and Dispensed Data
    order_pid <- data_order[data_order[,ID] == pid,]
    order_pid <- order_pid[!is.na(order_pid[,Drug]),]
    dd_pid <- data_dispensed[data_dispensed[,ID] == pid,]
    dd_pid <- dd_pid[!is.na(dd_pid[,Drug]),]

    # Use "ordering date" to substitute the missing "start date"
    for (i in 1:nrow(order_pid)) {
      if (is.na(order_pid[i,Start.Date])) {
        order_pid[i,Start.Date]=order_pid[i,Order.Date]
      }
    }

    # Use the date of data pulled / date of today to substitute the missing "end date"
    order_pid$end_date[is.na(order_pid$end_date)] <- Data.Date

    # Clinic Visit Data
    visit_pid <- data_visit[data_visit[,ID] == pid,]
    visits <- as.vector(as.matrix(unique(visit_pid[,Visit.Date])))
    death <- as.vector(as.matrix(unique(visit_pid[,Death.Date])))
    death <- as.Date(death, origin = "1970-01-01")

    df_pdc_id_visits <- NULL

    for (visit in visits) {

      visit <- as.Date(visit, origin = "1970-01-01")

      # Hospitalization Data
      if (sum(data_hosp[,ID] == pid) > 0){
        hosp_pid <- data_hosp[data_hosp[,ID] == pid,]
        hosp_pid <- unique(hosp_pid[,c(ID,Hosp.Start,Hosp.End)])
        # Use the date of data pulled to substitute the missing "discharge date"
        hosp_pid[,Hosp.End][is.na(Hosp.End),] <- Data.Date
      }
      if (sum(data_hosp[,ID] == pid) == 0){
        hosp_pid <- data.frame(matrix(ncol = 3))
        colnames(hosp_pid) <- c(ID,Hosp.Start,Hosp.End)
        hosp_pid[,ID] <- pid
      }

      ############################################################################
      # Find the denominator start date of each drug for the selected patient
      ############################################################################

      # Initializing variables and dataframe
      date_med_id <- NULL

      # Create data frame `df_pid_drugID` to summarize important info
      df_pid_drugID <- data.frame(matrix(ncol = 5))
      colnames(df_pid_drugID) <- c("ID",
                                   "Drug",
                                   "threshold_start",
                                   "threshold_end",
                                   "visit_date")

      # Create a new object `DrugIDs` for the drug list of the selected patient
      DrugIDs <- as.vector(as.matrix(unique(order_pid[,Drug])))
      drugID <- DrugIDs[1]

      ############################################################################
      # Obtain the start and end date of each medication
      ############################################################################

      for (drugID in DrugIDs) { # Loop over all DrugIDs
        # Create a new object `order_pid_drugID_cp` for orders of selected drug
        order_pid_drugID_cp <- as.data.frame(
          order_pid[order_pid[,Drug] == drugID,])

        # Sort the orders by start date and end date
        order_pid_drugID <- NULL
        order_pid_drugID$start_date <- sort(as.Date(
          order_pid_drugID_cp[,"start_date"]))
        order_pid_drugID$end_date <- sort(as.Date(
          order_pid_drugID_cp[,"end_date"]))
        order_pid_drugID <- as.data.frame(order_pid_drugID)

        ##########################################################################
        # Adjustment for Contiguous Medication Order:
        # the same medication was continued within a day of being discontinued
        ##########################################################################

        # Only one medication order: Interval NOT existed
        if (nrow(order_pid_drugID)==1) {
          df_pid_drugID$threshold_start <- order_pid_drugID$start_date
          df_pid_drugID$threshold_end <- order_pid_drugID$end_date}

        # Multiple medication orders: check the interval between orders
        if (nrow(order_pid_drugID)>1) {
          interval <- c()
          for (r in 1:(nrow(order_pid_drugID)-1)){
            interval[r] <- order_pid_drugID[,Start.Date][r+1]-
              order_pid_drugID[,End.Date][r]
          }
          interval[is.na(interval)]<-0

          if (all(interval<=1, na.rm=TRUE)) { # All orders are continual
            df_pid_drugID$threshold_start <- order_pid_drugID[,Start.Date][1]
            df_pid_drugID$threshold_end <- order_pid_drugID[,End.Date][nrow(
              order_pid_drugID)]
          }

          else{ # Interval existed between orders
            # Select the latest start date before the clinic visit date
            start_interval <- c(0,interval)
            end_interval <- c(interval,0)
            starts <- c(order_pid_drugID[,Start.Date][1])
            starts <- append(starts,
                             order_pid_drugID[,Start.Date][start_interval>1])
            # Create a logical vector `ind_order` to specify the order of date
            #ind_order = starts < visit
            #starts <- starts[ind_order]
            # Select the end date in line with the start date
            ends <- c(order_pid_drugID[,End.Date][end_interval>1])
            ends <- append(ends,
                           order_pid_drugID[,End.Date][nrow(order_pid_drugID)])
            #ends <- ends[ind_order]

            df_pid_drugID <- NULL
            for (i in 1:length(starts)) {
              mid <- data.frame(matrix(ncol = 5))
              colnames(mid) <- c("ID","Drug",
                                 "threshold_start","threshold_end","visit_date")
              mid$threshold_start <- starts[i]
              mid$threshold_end <- ends[i]
              df_pid_drugID <- rbind(df_pid_drugID, mid)
            }

            df_pid_drugID$threshold_start <- as.Date(df_pid_drugID$threshold_start)
            df_pid_drugID$threshold_end <- as.Date(df_pid_drugID$threshold_end)
          }
        }
        df_pid_drugID$ID <- order_pid_drugID_cp[,ID][1]
        df_pid_drugID$Drug <- order_pid_drugID_cp[,Drug][1]
        df_pid_drugID$visit_date <- visit

        # Return `date_med_id` as a dataframe with ID/Drug name/important Dates
        date_med_id <- rbind(date_med_id, df_pid_drugID)
      } # Close loop of all DrugIDs; Retune dataframe date_med_id

      date_med_id$threshold_start <- as.Date(date_med_id$threshold_start,
                                             origin="1970-01-01")
      date_med_id$threshold_end <- as.Date(date_med_id$threshold_end,
                                           origin="1970-01-01")
      date_med_id <- date_med_id %>%
        distinct() %>%
        filter(!is.na(threshold_start))

      # Create a column of logical constants to specify whether the medication is
      # active on the clinic visit date

      date_med_id <- date_med_id %>%
        #mutate(active_at_start = threshold_start <= Inclusion.Start & threshold_end > Inclusion.Start) %>%
        group_by(Drug) %>%
        mutate(first_prescription = min(threshold_start)) %>%
        mutate(end_prescription = max(threshold_end)) %>%
        mutate(if_new = first_prescription >= Inclusion.Start & first_prescription <= Inclusion.End) %>%
        mutate(exclude = first_prescription > Inclusion.End | end_prescription < Inclusion.Start) %>%
        filter(!exclude == TRUE)

      date_med_id <- as.data.frame(date_med_id)

      # Create new objects for the interested time frame of the selected patient
      if (all(date_med_id$first_prescription < Inclusion.Start)) {
        start <- Inclusion.Start
      }else {
        start <- max(date_med_id$first_prescription)
      }

      end <- start+Followed.Days
      total_days <- Followed.Days

      date_med_id <- date_med_id %>%
        mutate(start = start) %>%
        mutate(if_active = threshold_start<=start & threshold_end>start)


      if (!is.na(death) & death>=start & death<end) {
        end <- death
        total_days <- as.numeric(death-start)
      } # Loop closed for death

      # No active medications
      if (nrow(date_med_id) > 0) {

        ############################################################################
        # Find the dispensed dates for each medication order
        ############################################################################

        # Select dispensed data from 90 days prior 'start' to `end`
        dd_pid <- as.data.frame(dd_pid)
        dd_pid <- dd_pid[dd_pid[,Fill.Date] >= (start-90) &
                           dd_pid[,Fill.Date] < end,]

        dispense_pid <- NULL
        for (drugID in unique(date_med_id$Drug)) {
          dd_pid_drugID <- dd_pid[dd_pid[,Drug] == drugID,]

          # Sort dispensed data by fill date
          dd_pid_drugID <- dd_pid_drugID[order(dd_pid_drugID[,Fill.Date]),]
          dd_pid_drugID$Fill.End <- dd_pid_drugID[,Fill.Date] +
            dd_pid_drugID[,Days.Supplied]
          dispense_drugID <- dd_pid_drugID %>%
            dplyr::select(ID, Drug, Fill.Date, Fill.End) %>%
            rename(ID=ID) %>%
            rename(Drug=Drug) %>%
            rename(threshold_start=Fill.Date) %>%
            rename(threshold_end=Fill.End)

          dispense_pid <- rbind(dispense_pid,dispense_drugID)
        }
        # Return `dispense_pid` as a dataframe with start & end dates for each med

        ############################################################################
        # Find the stockpiled medications available
        ############################################################################

        # Create data frame `stockpiles` to summarize important info
        stockpiles <- data.frame(matrix(nrow = length(unique(date_med_id$Drug)),
                                        ncol = 3))
        row.names(stockpiles) <- unique(date_med_id$Drug)
        colnames(stockpiles) <- c("Stocked days", "Gap", "Gap_stockpile")
        date <- NULL
        stockpiles[,"Stocked days"] <- 0
        stockpiles[,"Gap"] <- 0
        stockpiles[,"Gap_stockpile"] <- 0

        for (dr in 1:length(unique(date_med_id$Drug))) {

          # Loop over all Drugs for the selected patient
          drugID <- unique(date_med_id$Drug)[dr]

          # Sort dispensed data by fill date
          dd_pid_drugID <- dd_pid[dd_pid[,Drug] == drugID,]
          dd_pid_drugID <- dd_pid_drugID[order(dd_pid_drugID[,Fill.Date]),]

          # Create a new object called `start_dr` for the selected medication
          starts_dr <- date_med_id[date_med_id[,"Drug"]==drugID &
                                     date_med_id[,"if_active"]==TRUE,]
          starts <- as.vector(as.matrix(unique(starts_dr[,"threshold_start"])))
          start_dr <- min(starts, na.rm = TRUE)
          start_dr <-  as.Date(start_dr, origin = "1970-01-01")
          start_dr <- max(start_dr, start)

          # Use logical index to get dispensed data for stockpiled

          ind_stockpiled <- (start_dr-90)<=dd_pid_drugID[,Fill.Date] &
            dd_pid_drugID[,Fill.Date]<start_dr

          ind_stockpiled[is.na(ind_stockpiled)] <- FALSE
          dates_stockpiled <- dd_pid_drugID[,Fill.Date][ind_stockpiled]
          first_stockpiled <- min(dates_stockpiled, na.rm = TRUE)
          period_stockpiled <- as.numeric(start-first_stockpiled)

          if (all(ind_stockpiled==FALSE,na.rm = TRUE)) {
            # NO dispensed data prior to the interested time frame
            days_stockpiled <- 0
          }
          else { # dispensed data existed prior to the interested time frame
            date_med_id %>%
              filter(Drug==drugID) -> od
            dispense_pid %>%
              filter(Drug==drugID) -> dd

            for (d in 0:(period_stockpiled-1)) {
              # Loop over the stockpiled time period
              date_s <- first_stockpiled + d
              # Use logical index to test whether the medication was ordered on the
              # specified date during the stockpiles period
              ind_date <- NULL
              for (i in 1:nrow(od)) {
                if (date_s >= od$threshold_start[i] &
                    (date_s <= od$threshold_end[i] | is.na(od$threshold_end[i]))) {
                  ind_d <- TRUE
                  ind_date <- append(ind_date, ind_d)
                }
                else {
                  ind_d <- FALSE
                  ind_date <- append(ind_date, ind_d)
                }
              }

              s <- 0
              gap_s <- stockpiles$Gap_stockpile[dr]

              ind_dispense <- NULL
              # Had order in this date (stockpiled)
              if (any(ind_date==TRUE)) {
                if (nrow(dd)>0) {
                  for (r in 1:nrow(dd)) {
                    if (date_s >= dd$threshold_start[r] & (date_s < dd$threshold_end[r] | is.na(dd$threshold_end[r]))) {
                      s <- s+1
                      ind_dd <- TRUE
                      ind_dispense <- append(ind_dispense,ind_dd)
                    }
                    else {
                      ind_dd <- FALSE
                      ind_dispense <- append(ind_dispense,ind_dd)}
                  }
                  # Had order but no dispense,stock-1
                }
                if (all(ind_dispense==FALSE) & stockpiles$`Stocked days`[dr]>0) {
                  stockpiles$`Stocked days`[dr] <- stockpiles$`Stocked days`[dr]-1
                }
                # Had not order in this date (stockpiled)
              }
              else {
                if (nrow(dd)>0) {
                  for (r in 1:nrow(dd)) {
                    if (date_s >= dd$threshold_start[r] & (date_s < dd$threshold_end[r] | is.na(dd$threshold_end[r]))) {
                      s <- 2
                    }
                  }
                }
                gap_s <- gap_s+1
              }
              if (s>1) {
                stockpiles$`Stocked days`[dr] <- stockpiles$`Stocked days`[dr]+(s-1)
              }
              if(gap_s > 30) {
                stockpiles$`Stocked days`[dr] <- 0
                gap_s <- 0
              }

              stockpiles$Gap_stockpile[dr] <- gap_s
              #print(date_s)
              #print(s)
              #print(stockpiles$`Stocked days`[dr])
            } # Close loop for the stockpiled time period
          }
        } # Close loop of all Drugs for the selected patient

        # stockpiles

        ##########################################################################
        # Covered days for each patient (numerator of PDC)
        ##########################################################################

        covered_days <- NULL
        covers <- data.frame(matrix(ncol = length(unique(date_med_id$Drug)),nrow = total_days))
        colnames(covers) <- c(unique(date_med_id$Drug))
        # sum_covers <- NULL
        dates <- NULL
        hosp_days <- NULL

        for (d in 0:(total_days-1)) {
          date <- start+d

          cover <- NULL

          for (dr in 1:length(unique(date_med_id$Drug))) {
            drug <- unique(date_med_id$Drug)[dr]
            dispense_pid %>%
              filter(Drug==drug) -> dd
            date_med_id %>%
              filter(Drug==drug) -> od

            #
            ind_date <- NULL
            for (i in 1:nrow(od)) {
              if (date >= od$threshold_start[i] & (date <= od$threshold_end[i] | is.na(od$threshold_end[i]))) {
                ind_d <- TRUE
                ind_date <- append(ind_date,ind_d)
              }
              else {
                ind_d <- FALSE
                ind_date <- append(ind_date,ind_d)
              }
            }

            # Check if this date has order for this med, return ind_date


            ############################################################
            # Hospitalization

            hosp_pid <- as.data.frame(hosp_pid)
            ind_hosp <- NULL
            if (nrow (hosp_pid)>1) {
              for (i in 1:nrow(hosp_pid)) {
                if (date >= hosp_pid[,Hosp.Start][i] & (date <= hosp_pid[,Hosp.End][i] | is.na(hosp_pid[,Hosp.End][i])))
                {
                  ind_h <- TRUE
                  ind_hosp <- append(ind_hosp,ind_h)
                }
                else {
                  ind_h <- FALSE
                  ind_hosp <- append(ind_hosp,ind_h)
                }
              }
            }
            # Return ind_hosp for if hospitalized on this date

            if (any(ind_hosp==TRUE)) {
              in_hosp <- TRUE
              ind_date = FALSE
            } else {
              in_hosp <- FALSE
            }
            hosp_days <- append(hosp_days,in_hosp)

            ############################################################
            gap <- stockpiles$Gap[dr]
            if (any(ind_date==TRUE)) {
              c <- 0
              gap <- 0
              if (nrow(dd)>0) {
                for (r in 1:nrow(dd)) {
                  if (date >= dd$threshold_start[r] & (date < dd$threshold_end[r] | is.na(dd$threshold_end[r]))) {
                    c <- c+1
                  }
                }
              }
            }else {c <- NA
            gap <- gap+1
            }

            adjust <- stockpiles$`Stocked days`[dr]
            if (is.na(c)) {
              if (gap > 30) {
                adjust <- 0
              }
              if (gap <= 30) {
                NULL
              }
            }else if (c>1) {
              adjust <- adjust+1
            }else if(c==0 & adjust>0) {
              c <- c+1
              adjust <- adjust-1
            }

            stockpiles$`Stocked days`[dr] <- adjust
            stockpiles$Gap[dr] <- gap

            ### cover <- append(cover,c)
            #print(paste(drug,": ",date,": ",gap,": ",adjust))
            #sum <- sum(!is.na(cover)) == sum(cover>0,na.rm = TRUE)
            #if(all(is.na(cover))) {sum <- FALSE}
            #sum_covers <- append(cover,sum)
            #covers[,d+1] <- sum_covers
            ### !!! return c as number of pills of this med on that date

            cover <- append(cover,c)
            #print(dr)
            #print(cover)
            #print(dates)
          }

          dates <- append(dates,date)
          covers[d+1,] <- cover
          total_hosp_days <- sum(hosp_days)
        }

        row.names(covers) <- dates
        # Return covers as a table of all meds
        # covered_days <- sum(covers["Sum",]==TRUE)

        total_hosp_days ### Total hospitalization days in these 180 days


        ###################################
        # Calculate PDC
        ###################################

        # PDC <- NULL
        # To remove any gap in denominator
        #gap_count <- covers[1:length(unique(date_med_id$Drug)),]
        #gap_count[,'Gap_count'] <- apply(gap_count, MARGIN = 1, function(x) sum(is.na(x))) == length(unique(date_med_id$Drug))
        #gap_denomintor <- sum(gap_count['Gap_count',]==TRUE)
        # Calculating PDC

        # total_days <- total_days - gap_denomintor
        # PDC <- as.numeric(covered_days)/total_days

        # Set PDC to 1 if PDC is over 1
        #if(PDC>1 & !is.na(PDC)) {df_pid_drugID$pdc=1}

        # Set minimum threshold of 30 days for denominator
        # df_pdc_id <- subset(df_pdc_id,total_days>=30)

        ##########################################################################
        # Results
        ##########################################################################

        df_pdc_id <- NULL

        df_pdc_id_visit <- data.frame(matrix(ncol = 4))
        colnames(df_pdc_id_visit) <- c("ID","threshold_start","threshold_end","visit_date")
        df_pdc_id_visit$ID <- pid
        df_pdc_id_visit$threshold_start <- start
        df_pdc_id_visit$threshold_end <- end
        df_pdc_id_visit$visit_date <- visit
        df_pdc_id_visit$death_date <- death
        df_pdc_id_visit$hospitalization_days <- total_hosp_days


        covers <- covers %>%
          rowwise() %>%
          mutate(active_meds = sum(!is.na(c_across(1:length(unique(date_med_id$Drug)))))) %>%
          mutate(on_meds = sum(c_across(1:length(unique(date_med_id$Drug))),na.rm = TRUE)) %>%
          mutate(minimum = ifelse(active_meds>0 & on_meds==active_meds,TRUE,ifelse(active_meds==0,NA,FALSE))) %>%
          mutate(maximum = ifelse(active_meds>0 & on_meds>=1,TRUE,ifelse(active_meds==0,NA,FALSE)))


        min_covered <- sum(covers[,"minimum"]==TRUE, na.rm = TRUE)
        max_covered <- sum(covers[,"maximum"]==TRUE, na.rm = TRUE)

        total <- sum(covers[,"active_meds"]>0)
        df_pdc_id_visit$total_days <- total

        df_pdc_id_visit$minimum_pdc <- min_covered/total
        df_pdc_id_visit$maximum_pdc <- max_covered/total

        df_pdc_id_visit$gap <- sum(covers[,"active_meds"]==0)


        # Need to remove it if denominator=0??
        # df_pdc_id <- df_pdc_id %>%
        # filter(!is.na(pdc))

        df_pdc_id <- rbind(df_pdc_id,df_pdc_id_visit)
      } else { # Loop closed for no active meds
        df_pdc_id <- data.frame(matrix(ncol = 4))
        colnames(df_pdc_id) <- c("ID","threshold_start","threshold_end","visit_date")
        df_pdc_id$ID <- pid
        df_pdc_id$threshold_start <- start
        df_pdc_id$threshold_end <- end
        df_pdc_id$visit_date <- visit
        df_pdc_id$death_date <- death
        df_pdc_id$hospitalization_days <- NA
        df_pdc_id$total_days <- NA
        df_pdc_id$minimum_pdc <- NA
        df_pdc_id$maximum_pdc <- NA
        df_pdc_id$gap <- NA } # Loop closed for no active meds else

      # Return df_pdc_id

      df_pdc_id_visits <- rbind(df_pdc_id_visits,df_pdc_id)

    }# Loop closed for visits

    return(df_pdc_id_visits)
  }

  ##############################################################################
  # End of the inner function
  ##############################################################################

  pids <- as.vector(as.matrix(unique(data_order[,ID])))
  patientID <- 1:length(pids)
  df_pdc <- NULL

  # Whether to use parallel to calculate for all patients
  if(parallel == T){

    clnum <- detectCores()
    if( clnum <=  ncores ){
      clnum <- clnum
      print(paste('only', clnum, 'cores was detected, change the number of cores to', clnum))

      #Sys.sleep(1)
    }else{
      clnum < ncores
    }

    #clnum <- 4
    cl <- makeCluster(clnum)
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = length(patientID), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    df_pdc <- foreach(x=patientID,.packages=c('lubridate', 'tidyverse'),
                      .options.snow = opts,
                      .combine = rbind ) %dopar% {

                     pdc_min_max_pid(data_dispensed = data_dispensed,
                                     patientID = x,
                                     ID = ID,
                                     Drug = Drug,
                                     Fill.Date = Fill.Date,
                                     Days.Supplied = Days.Supplied,
                                     data_order = data_order,
                                     Order.Date=Order.Date,
                                     Start.Date=Start.Date,
                                     End.Date=End.Date,
                                     data_visit = data_visit,
                                     Visit.Date = Visit.Date,
                                     Death.Date = Death.Date,
                                     Followed.Days = 180,
                                     data_hosp = data_hosp,
                                     Hosp.Start = Hosp.Start,
                                     Hosp.End = Hosp.End,
                                     Inclusion.Start = Inclusion.Start,
                                     Inclusion.End = Inclusion.End)
                      }

  }else{
    for (pid in patientID) {
      pdc_id <- pdc_min_max_pid(data_dispensed = data_dispensed,
                             patientID = pid,
                             ID = ID,
                             Drug = Drug,
                             Fill.Date = Fill.Date,
                             Days.Supplied = Days.Supplied,
                             data_order = data_order,
                             Order.Date = Order.Date,
                             Start.Date = Start.Date,
                             End.Date = End.Date,
                             data_visit = data_visit,
                             Visit.Date = Visit.Date,
                             Death.Date = Death.Date,
                             Followed.Days = 180,
                             data_hosp = data_hosp,
                             Hosp.Start = Hosp.Start,
                             Hosp.End = Hosp.End,
                             Inclusion.Start = Inclusion.Start,
                             Inclusion.End = Inclusion.End)

      df_pdc <- rbind(df_pdc, pdc_id)
      cat('\014')
      print(paste0(round(pid / length(patientID),4) * 100, '%'))
    }

    df_pdc
  }

  return(df_pdc)
}
