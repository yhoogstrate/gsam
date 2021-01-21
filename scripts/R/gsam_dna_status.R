#!/usr/bin/env R 

tmp <- read.csv('data/administratie/GSAM_combined_clinical_molecular.csv', stringsAsFactors=F)
rownames(tmp) <- tmp$studyID
tmp$X.1 <- NULL
tmp$X <- NULL
#tmp$studyID <- NULL

# not needed for this sudy?
tmp$gender <- NULL
tmp$age <- NULL
tmp$status <- NULL
tmp$survivalDays <- NULL
tmp$nDocumentedRecurrences <- NULL
tmp$progressionFreeDays <- NULL
tmp$survivalFromSecondSurgeryDays <- NULL
tmp$extentOfResectionFirstSurgery <- NULL
tmp$extentFirstSurgeryConfirmedByMRI <- NULL
tmp$tumorLocation <- NULL
tmp$initialHistology <- NULL
tmp$performanceAtStartRT <- NULL
tmp$steroidsAtStartRT <- NULL
tmp$steroidsDosageAtStartRT <- NULL
tmp$initialMGMT <- NULL
tmp$treatedWithRT <- NULL
tmp$dosageRT <- NULL
tmp$treatedWithTMZ <- NULL
tmp$daysToCompletionRT <- NULL
tmp$nCyclesTMZ <- NULL
tmp$daysToLastCycleTMZ <- NULL
tmp$reasonDiscontinuationTMZ <- NULL
tmp$otherTreatmentsBeforeSecondSurgery <- NULL
tmp$daysToSecondSurgery <- NULL
tmp$secondSurgeryForRecurrenceNumber <- NULL
tmp$extentOfResectionSecondSurgery <- NULL
tmp$extentSecondSurgeryConfirmedByMRI <- NULL
tmp$recurrentHistology <- NULL
tmp$performanceAtSecondSurgery <- NULL
tmp$sizeLargestLesionSecondSurgery <- NULL
tmp$nTargetLesionsSecondSurgery <- NULL
tmp$steroidsBeforeSecondSurgery <- NULL
tmp$treatmentDetailsFirstPD <- NULL
tmp$discontinuationReasonTreatmentFirstPD <- NULL
tmp$daysToSecondPD <- NULL
tmp$treatmentForSecondPD <- NULL
tmp$treatmentDetailsSecondPD <- NULL
tmp$discontinuationReasonTreatmentSecondPD <- NULL
tmp$daysToThirdPD <- NULL
tmp$treatmentForThirdPD <- NULL
tmp$treatmentDetailsThirdPD <- NULL
tmp$discontinuationReasonTreatmentThirdPD <- NULL
tmp$daysToFourthPD <- NULL
tmp$treatmentForFourthPD <- NULL
tmp$treatmentDetailsFourthPD <- NULL
tmp$discontinuationReasonTreatmentFourthPD <- NULL
tmp$daysToFifthPD <- NULL
tmp$centerLetter <- NULL
tmp$PID <- NULL
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





