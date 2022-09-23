#!/usr/bin/env R


# useful link: https://emilhvitfeldt.github.io/r-color-palettes/discrete.html

mixcol <- function(c1, c2, ratio=0.5) {
  n <- names(c1)

  m1 <- col2rgb(c1)
  m2 <- col2rgb(c2)
  
  m <- (ratio * m2) + ((1 - ratio) * m1)
  r <- rgb(m[1,],m[2,],m[3,], maxColorValue = 255)
  names(r) <- n

  return(r)
}

# test func
#mixcol(subtype_colors_nmf2022, rep("black",length(subtype_colors_nmf2022)),0.15)



subtype_colors <- c('Classical'='#6ba6e5',# blue
                    'Mesenchymal'='#eab509',#mustard
                    'Proneural'='#ff5f68')#red/pink

subtype_colors_ext <- subtype_colors
subtype_colors_ext['Proneural|Classical'] <- mixcol(subtype_colors_ext['Proneural'] , subtype_colors_ext['Classical'] )
#subtype_colors_ext['Classical|Proneural'] <- mixcol(subtype_colors_ext['Proneural'] , subtype_colors_ext['Classical'] )

subtype_colors_nmf2022 <- c('NMF meta-feature 1'='#6ba6e5',# blue
                            'NMF meta-feature 3'='#eab509',#mustard
                            'NMF meta-feature 2'='#ff5f68')#red/pink

subtype_colors_ssGSEA <- c('ssGSEA CL score'='#6ba6e5',# blue
                            'ssGSEA MES score'='#eab509',#mustard
                            'ssGSEA PN score'='#ff5f68')#red/pink



subtype_colors_nmf <- c('NMF:123456.2'='#6ba6e5',# blue
                        'NMF:123456.1'='#eab509',#mustard
                        'NMF:123456.3'='#ff5f68')#red/pink


subtype_colors_nmfm <- c('NMF meta-feature 2'='#6ba6e5',# blue
                         'NMF meta-feature 1'='#eab509',#mustard
                         'NMF meta-feature 3'='#ff5f68')#red/pink

dataset_colors <- c('G-SAM' = '#e69356', 'GSAM' = '#e69356', 'GLASS' = '#69a4d5')

resection_colors <- c(
  "primary" = "#55864e",
  "Primary Res." = "#55864e",
  "R1" = "#55864e" ,
  
  "recurrence" = "#a15e5e",
  "Recurrent Res." = "#a15e5e",
  "R2" = "#a15e5e"  # #a1805e
)

