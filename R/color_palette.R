color_pal <- function(n = 20, color, trim_upp, trim_low){
  pal <- hcl.colors(20, "Blues")[-c(1:2)]
  show_image(pal)
  as.matrix(t(col2rgb(pal)))
}
