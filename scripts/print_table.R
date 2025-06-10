print_table <- function(df,
                        title = NULL,
                        row.names = F,
                        font_family = "serif",
                        save.pdf = F,
                        path2save,
                        file.name = "table",
                        width = 10,
                        height = 10) {
  library(gridExtra)
  library(kableExtra)
  library(gtable)
  
  print(
  df %>%
    kbl(caption = title, row.names = row.names) %>%
    kable_classic(full_width = F, html_font = font_family)
  )
  
  if (save.pdf) {
    pdf(paste0(path2save,"/",file.name,".pdf"), width = width, height = height)
    tt = ttheme_default(base_family = font_family,
                        core=list(
                          bg_params = list(fill=c(rep(c("grey95", "grey90"),
                                                      length.out=nrow(df)))
                          )))
    grid.table(df, theme = tt, rows = rep("",nrow(df)))
    dev.off()
  }
  
}