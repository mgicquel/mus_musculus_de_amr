
# this script contains the premade themes for the ggplots for three options:
# - theme_mg: light background single plot
# - theme_mg_facets: light background facet plot
# - theme_mg_dark: dark background single plot


# to be able to change the fonts, we might need to run this code before
library(extrafont)
#font_import()
loadfonts(device = "win")

theme_mg <- function(){

  font <- "Poppins"   #assign font family up front

  theme_minimal() %+replace%    #replace elements we want to change

    theme(

      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_line(colour = "black"),          #strip axis ticks
      axis.line = element_line(colour = "black"),           # axis lines

      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 10,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly

      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 9),               #font size

      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align

      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 9,
        face = "bold"),               #font size

      axis.text = element_text(              #axis text
        family = font,            #axis family
        size = 8),                #font size

      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10)),

      axis.text.y = element_text(            #margin for axis text
        margin=margin(l = 5, r = 10)),

      legend.title = element_text(             #axis titles
        family = font,            #font family
        size = 10,
        face = "bold"),

      legend.text = element_text(              #axis text
        family = font,            #axis family
        size = 9)

      #since the legend often requires manual tweaking
      #based on plot content, don't define it here
    )
}



theme_mg_facets <- function(){

  font <- "Poppins"  #assign font family up front

  theme_minimal() %+replace%    #replace elements we want to change

    theme(

      #grid elements
      # panel.grid.major = element_blank(),    #strip major gridlines
      # panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_line(colour = "black"),          #strip axis ticks
      axis.line = element_line(colour = "black"),           # axis lines

      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 12,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly

      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 12),               #font size

      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align

      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 12,
        face = "bold"),               #font size

      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 10),                #font size

      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10)),

      axis.text.y = element_text(            #margin for axis text
        margin=margin(l = 5, r = 10)),

      legend.title = element_text(             #axis titles
        family = font,            #font family
        size = 10,
        face = "bold"),

      legend.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 10),

      strip.text.x = element_text(
        family = font,
        size = 10, color = "black"),

      strip.text.y = element_text(
        family = font,
        size = 10, color = "black")

    )
}


### THEME FOR MAP

theme_mg_map <- function(){

  font <- "Poppins"   #assign font family up front

  theme_bw() %+replace%    #replace elements we want to change

    theme(


      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_line(colour = "black"),          #strip axis ticks
      axis.line = element_line(colour = "black"),           # axis lines

      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 12,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly

      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 10),               #font size

      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align

      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 11,
        face = "bold"),               #font size

      axis.text = element_text(              #axis text
        family = font,            #axis family
        size = 9),                #font size

      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10)),

      axis.text.y = element_text(            #margin for axis text
        margin=margin(l = 5, r = 10)),

      legend.title = element_text(             #axis titles
        family = font,            #font family
        size = 11,
        face = "bold"),

      legend.text = element_text(              #axis text
        family = font,            #axis family
        size = 9)

      #since the legend often requires manual tweaking
      #based on plot content, don't define it here
    )
}


### THEME FOR SANKEY DIAGRAM

theme_mg_sankey <- function(){

  font <- "Poppins"  #assign font family up front

  theme_sankey(base_family = "Poppins") %+replace%    #replace elements we want to change

    theme(

      #grid elements
      # panel.grid.major = element_blank(),    #strip major gridlines
      # panel.grid.minor = element_blank(),    #strip minor gridlines
      #axis.ticks = element_line(colour = "black"),          #strip axis ticks
      #axis.line = element_line(colour = "black"),           # axis lines

      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 20,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly

      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 14),               #font size

      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align

      axis.title = element_blank(),               #font size

      axis.text = element_text(              #axis text
        family = font,            #axis family
        size = 10),                #font size

      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10)),

      axis.text.y = element_blank(),

      legend.title = element_text(             #axis titles
        family = font,            #font family
        size = 12,
        face = "bold"),

      legend.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 10)

    )
}


### THEME FOR CORRELATION PLOT

theme_mg_cor <- function(){

  font <- "Poppins"   #assign font family up front

  theme_minimal() %+replace%    #replace elements we want to change

    theme(

      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_blank(),          #strip axis ticks
      axis.line = element_blank(),           # axis lines

      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 12,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly

      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 10),               #font size

      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align

      axis.title = element_blank(),               #font size

      axis.text = element_text(              #axis text
        family = font,            #axis family
        size = 9),                #font size

      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10), angle = 45, hjust = 1, vjust = 1),

      axis.text.y = element_text(            #margin for axis text
        margin=margin(l = 5, r = 10)),

      legend.title = element_text(             #axis titles
        family = font,            #font family
        size = 11,
        face = "bold"),

      legend.text = element_text(              #axis text
        family = font,            #axis family
        size = 9)

      #since the legend often requires manual tweaking
      #based on plot content, don't define it here
    )
}

### THEME FOR CIRCULAR CORRELATION PLOT

theme_mg_cor_cir <- function(){

  font <- "Poppins"   #assign font family up front

  theme_minimal() %+replace%    #replace elements we want to change

    theme(

      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_blank(),          #strip axis ticks
      axis.line = element_blank(),           # axis lines
      plot.background = element_blank(),   # No plot background
      panel.background = element_blank(),  # No panel background

      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 12,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly

      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 10),               #font size

      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align

      axis.title = element_blank(),               #font size

      axis.text = element_blank(),                #font size

      legend.title = element_text(             #axis titles
        family = font,            #font family
        size = 11,
        face = "bold"),

      legend.text = element_text(              #axis text
        family = font,            #axis family
        size = 9),

      # Positioning and orientation of the legend
      legend.position = "bottom",          # Legend at the bottom
      legend.direction = "horizontal",     # Horizontal legend
      legend.box = "vertical",             # Vertically stacked guides

      #since the legend often requires manual tweaking
      #based on plot content, don't define it here
    )
}

