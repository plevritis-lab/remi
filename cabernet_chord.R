
# Chord Diagram function
countLRs <- function(mat, combo) {
  count <- unique(pull(mat %>% filter(celltypes == combo) %>%
                         select(count)))
  if(length(count) == 0) {
    return(0)
  } else {
    return(count)
  }
}

# Code to draw chord diagram. I actually manually put in the numbers (pretty inefficient) for the output 
# in illustrator because I wanted to control the position
# Sig.maxedges = output of the cabernet function
png("NAME_OF_IMAGE.png", width = 6, height = 6, units="in", res=600)
chord.sigmaxedges <- sig.maxedges %>%
  mutate(node1 = L) %>%
  mutate(node2 = R) %>%
  dplyr::select(node1, node2, n1cell, n2cell) %>%
  unique()

strsplit.ind <- seq(from=1, by=2, length.out = nrow(chord.sigmaxedges))
adeno.lr <- chord.sigmaxedges %>% 
  mutate(celltypes = paste0(n1cell, "_", n2cell)) %>%
  group_by(celltypes) %>%
  mutate(count = n()) %>%
  ungroup() 

# Sorry this is hardcoded! It's just all the combinations of the cell types which you can probably write a for loop to do
df <- data.frame(from=c('Fibroblast', 'Immune', 'Endothelial', 
                        'Immune', 'Cancer', 'Endothelial',
                        'Endothelial', 'Cancer', 'Fibroblast', 
                        'Fibroblast', 'Immune', 'Cancer',
                        'Fibroblast', 'Immune', 'Endothelial', 'Cancer'),
                 to=c('Cancer', 'Cancer', 'Cancer', 
                      'Fibroblast', 'Fibroblast', 'Fibroblast',
                      'Immune', 'Immune', 'Immune', 
                      'Endothelial', 'Endothelial', 'Endothelial',
                      'Fibroblast', 'Immune', 'Endothelial', 'Cancer'),
                 value=c(countLRs(adeno.lr, "F_M"),
                         countLRs(adeno.lr, "I_M"),
                         countLRs(adeno.lr, "E_M"),
                         countLRs(adeno.lr, "I_F"),
                         countLRs(adeno.lr, "M_F"),
                         countLRs(adeno.lr, "E_F"),
                         countLRs(adeno.lr, "E_I"),
                         countLRs(adeno.lr, "M_I"),
                         countLRs(adeno.lr, "F_I"),
                         countLRs(adeno.lr, "F_E"),
                         countLRs(adeno.lr, "I_E"),
                         countLRs(adeno.lr, "M_E"),
                         countLRs(adeno.lr, "F_F"),
                         countLRs(adeno.lr, "I_I"),
                         countLRs(adeno.lr, "E_E"),
                         countLRs(adeno.lr, "M_M")), stringsAsFactors=FALSE)
# Colors for the different cell types
grid.col <- c(Fibroblast='#DFC27D', 
              Cancer="#44AA99", 
              Immune="#99CCEE", 
              Endothelial="#A06846")
lwd_mat = matrix(1, nrow = nrow(df), ncol = ncol(df))
# Comment out code below allows you to control the thickness of the boundaries
# in the chord diagram
#lwd_mat[df$from == 'Cancer'] = 10
chordDiagram(df, grid.col=grid.col, 
             directional=1,direction.type = c("diffHeight", "arrows"), 
             annotationTrack = c("grid"), link.lwd = lwd_mat) 
dev.off()
