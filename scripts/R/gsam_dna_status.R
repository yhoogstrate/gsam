#!/usr/bin/env R 

# load libs ----


library(tidyverse)


# load data ----


#tmp$studyID <- NULL
tmp <- read.csv('data/gsam/administratie/GSAM_combined_clinical_molecular.csv', stringsAsFactors=F) %>%
  tibble::column_to_rownames('studyID') %>%
  dplyr::mutate(X.1 = NULL, X = NULL) %>%
  dplyr::mutate(gender = NULL) %>%
  dplyr::mutate(age = NULL) %>%
  dplyr::mutate(status = NULL) %>%
  dplyr::mutate(survivalDays = NULL) %>%
  dplyr::mutate(nDocumentedRecurrences = NULL) %>%
  dplyr::mutate(progressionFreeDays = NULL) %>%
  dplyr::mutate(survivalFromSecondSurgeryDays = NULL) %>%
  dplyr::mutate(extentOfResectionFirstSurgery = NULL) %>%
  dplyr::mutate(extentFirstSurgeryConfirmedByMRI = NULL) %>%
  dplyr::mutate(tumorLocation = NULL) %>%
  dplyr::mutate(initialHistology = NULL) %>%
  dplyr::mutate(performanceAtStartRT = NULL) %>%
  dplyr::mutate(steroidsAtStartRT = NULL) %>%
  dplyr::mutate(steroidsDosageAtStartRT = NULL) %>%
  dplyr::mutate(initialMGMT = NULL) %>%
  dplyr::mutate(treatedWithRT = NULL) %>%
  dplyr::mutate(dosageRT = NULL) %>%
  dplyr::mutate(treatedWithTMZ = NULL) %>%
  dplyr::mutate(daysToCompletionRT = NULL) %>%
  dplyr::mutate(nCyclesTMZ = NULL) %>%
  dplyr::mutate(daysToLastCycleTMZ = NULL) %>%
  dplyr::mutate(reasonDiscontinuationTMZ = NULL) %>%
  dplyr::mutate(otherTreatmentsBeforeSecondSurgery = NULL) %>%
  dplyr::mutate(daysToSecondSurgery = NULL) %>%
  dplyr::mutate(secondSurgeryForRecurrenceNumber = NULL) %>%
  dplyr::mutate(extentOfResectionSecondSurgery = NULL) %>%
  dplyr::mutate(extentSecondSurgeryConfirmedByMRI = NULL) %>%
  dplyr::mutate(recurrentHistology = NULL) %>%
  dplyr::mutate(performanceAtSecondSurgery = NULL) %>%
  dplyr::mutate(sizeLargestLesionSecondSurgery = NULL) %>%
  dplyr::mutate(nTargetLesionsSecondSurgery = NULL) %>%
  dplyr::mutate(steroidsBeforeSecondSurgery = NULL) %>%
  dplyr::mutate(treatmentDetailsFirstPD = NULL) %>%
  dplyr::mutate(discontinuationReasonTreatmentFirstPD = NULL) %>%
  dplyr::mutate(daysToSecondPD = NULL) %>%
  dplyr::mutate(treatmentForSecondPD = NULL) %>%
  dplyr::mutate(treatmentDetailsSecondPD = NULL) %>%
  dplyr::mutate(discontinuationReasonTreatmentSecondPD = NULL) %>%
  dplyr::mutate(daysToThirdPD = NULL) %>%
  dplyr::mutate(treatmentForThirdPD = NULL) %>%
  dplyr::mutate(treatmentDetailsThirdPD = NULL) %>%
  dplyr::mutate(discontinuationReasonTreatmentThirdPD = NULL) %>%
  dplyr::mutate(daysToFourthPD = NULL) %>%
  dplyr::mutate(treatmentForFourthPD = NULL) %>%
  dplyr::mutate(treatmentDetailsFourthPD = NULL) %>%
  dplyr::mutate(discontinuationReasonTreatmentFourthPD = NULL) %>%
  dplyr::mutate(daysToFifthPD = NULL) %>%
  dplyr::mutate(centerLetter = NULL) %>%
  dplyr::mutate(PID = NULL)


tmp[is.na(tmp)] <- "N/A"



# decompose from paired to single
tmp.p <- tmp
rownames(tmp.p) <- paste0(rownames(tmp.p), '1')
tmp.p[tmp.p == "Stable"] <- "Mutant"
tmp.p[tmp.p == "Lost"] <- "Mutant"
tmp.p[tmp.p == "Gained"] <- "Wildtype"


tmp.r <- tmp
rownames(tmp.r) <- paste0(rownames(tmp.r), '2')
tmp.r[tmp.r == "Stable"] <- "Mutant"
tmp.r[tmp.r == "Lost"] <- "Wildtype"
tmp.r[tmp.r == "Gained"] <- "Mutant"


gsam.dna.status <- rbind(tmp.p, tmp.r)
gsam.dna.status[gsam.dna.status == "N/A"] <- NA
gsam.dna.status$sid <- rownames(gsam.dna.status)

rm(tmp, tmp.p, tmp.r)





