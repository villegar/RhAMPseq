hex_logo <- function(subplot = "figures/slope.png",
                     dpi = 800,
                     h_color = "#000000",
                     h_fill = "#363b74",
                     output = "figures/metapipe.png",
                     p_color = "#eeeeee") {
  hexSticker::sticker(subplot = subplot, package = "RhAMPseq", 
                      h_color = h_color,  h_fill = h_fill,
                      dpi = dpi,
                      # l_x = 1.0, l_y = 1.0, spotlight = FALSE, 
                      s_x = 1.0, s_y = .85, s_width = .5,
                      p_x = 1.0, p_y = 1.52, p_size = 6, p_color = p_color,
                      url = "https://github.com/villegar/RhAMPseq", 
                      u_angle = 30, u_color = p_color, u_size = 1.35,
                      filename = output)
}

hex_logo(dpi = 500, h_fill = "#4F0924", output = "figures/logo.png")
