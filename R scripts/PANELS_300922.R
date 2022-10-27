library(ggpubr)
library(ggplot2)
library(cowplot)

#FS12
ggarrange(FS12A,FS12B,FS12C,ncol = 1, nrow = 3, common.legend = T)

#FS11
ggarrange(annotate_figure(FS11A, fig.lab = 'a', fig.lab.size = 18),
          annotate_figure(FS11B, fig.lab = 'b', fig.lab.size = 18),
          annotate_figure(FS11C, fig.lab = 'c', fig.lab.size = 18),
          ncol = 1, nrow = 3, heights = c(6,2,2.3))

#F2
ggarrange(annotate_figure(F2A, fig.lab = 'a', fig.lab.size = 18),
          annotate_figure(F2B, fig.lab = 'b', fig.lab.size = 18),
          annotate_figure(F2C, fig.lab = 'c', fig.lab.size = 18),
          annotate_figure(F2D, fig.lab = 'd', fig.lab.size = 18),
          ncol = 2, nrow = 2, widths = c(3,2.5))

#F6
ggarrange(annotate_figure(F6A, fig.lab = 'a', fig.lab.size = 18),
          annotate_figure(plot_grid(F6B, F6C, align = 'h', axis = 'tb'), fig.lab = 'b', fig.lab.size = 18),
          ggarrange(annotate_figure(F6D, fig.lab = 'd', fig.lab.size = 18),
                    annotate_figure(F6E, fig.lab = 'e', fig.lab.size = 18), ncol = 2, nrow = 1),
          ncol = 1, nrow = 3, heights = c(1.25,1.5,1.5))

#FS6
ggarrange(annotate_figure(FS6A, fig.lab = 'a', fig.lab.size = 18),
          annotate_figure(FS6B, fig.lab = 'b', fig.lab.size = 18),
          annotate_figure(FS6C, fig.lab = 'c', fig.lab.size = 18),
          ncol = 3, nrow = 1, widths = c(2,2,0.75))

#FS14
ggarrange(
ggarrange(annotate_figure(FS14A, fig.lab = 'a', fig.lab.size = 18),
          annotate_figure(FS14B, fig.lab = 'b', fig.lab.size = 18), nrow = 1, widths = c(1.75,1)),
          annotate_figure(FS14C, fig.lab = 'c', fig.lab.size = 18),
          ncol = 1, nrow = 2, heights = c(1,2))

#F8
ggarrange(
ggarrange(annotate_figure(F8A, fig.lab = 'a', fig.lab.size = 18),
          annotate_figure(F8B, fig.lab = 'b', fig.lab.size = 18), ncol = 1, nrow = 2, heights = c(1,0.75)),
          annotate_figure(F8C, fig.lab = 'c', fig.lab.size = 18),
  ncol = 2, nrow = 1, widths = c(1.1,1))

#FS2
ggarrange(annotate_figure(FS2B, fig.lab = 'b', fig.lab.size = 18),
ggarrange(annotate_figure(FS2C, fig.lab = 'c', fig.lab.size = 18),
          annotate_figure(FS2D, fig.lab = 'd', fig.lab.size = 18), nrow=1, widths = c(0.46,0.54)),
          annotate_figure(FS2E, fig.lab = 'e', fig.lab.size = 18), nrow = 3, heights = c(0.75,1,0.75))

#FS7
ggarrange(
  annotate_figure(plot_grid(FS7A, FS7C, align = 'h', axis = 'lr'), fig.lab.size = 18, fig.lab = 'a'),
  ggarrange(annotate_figure(FS7B, fig.lab.size = 18, fig.lab = 'b'), 
            annotate_figure(FS7D, fig.lab.size = 18, fig.lab = 'd'), nrow = 1),
  annotate_figure(plot_grid(FS7E, FS7G, align = 'h', axis = 'lr'), fig.lab.size = 18, fig.lab = 'e'),
  ggarrange(annotate_figure(FS7F, fig.lab.size = 18, fig.lab = 'f'), 
            annotate_figure(FS7H, fig.lab.size = 18, fig.lab = 'h'), nrow = 1), 
  ncol = 1, heights = c(0.7,0.3,1,0.3))

#F4
ggarrange(
  annotate_figure(plot_grid(F4A, F4B, F4C, align = 'h', axis = 'lr', nrow = 1, rel_widths = c(5, 5, 20)),
                fig.lab = 'abc', fig.lab.size = 18),
  annotate_figure(F4D, fig.lab = 'd', fig.lab.size = 18), 
nrow = 2, heights = c(0.75,1))

#FS9
ggarrange(annotate_figure(FS9A, fig.lab = 'a', fig.lab.size = 18),
          annotate_figure(FS9B, fig.lab = 'b', fig.lab.size = 18),
          nrow = 1)

#F5
ggarrange(
ggarrange(annotate_figure(F5A, fig.lab = 'a', fig.lab.size = 18),
          annotate_figure(F5B, fig.lab = 'b', fig.lab.size = 18),
          nrow = 1, widths = c(0.75,0.3)),
annotate_figure(F5C, fig.lab = 'c', fig.lab.size = 18), nrow = 2, heights = c(0.3,0.65))

#F3
ggarrange(annotate_figure(F3A+theme(plot.margin = margin(l=60,r=60)), fig.lab = 'a', fig.lab.size = 18),
          annotate_figure(F3B+theme(plot.margin = margin(l=60,r=60)), fig.lab = 'b', fig.lab.size = 18),
          annotate_figure(F3C, fig.lab = 'c', fig.lab.size = 18),
          nrow = 3, heights = c(0.7,0.3,1))

#F1
ggarrange(annotate_figure(F1A, fig.lab.size = 18, fig.lab = 'a'),
          annotate_figure(plot_grid(F1B,F1C,F1D,align = 'v', axis = 'lr', nrow = 3), 
                          fig.lab.size = 18, fig.lab = 'bcd', fig.lab.pos = 'top'),
          nrow = 1, widths = c(0.4,0.6))

#F7
ggarrange(
ggarrange(annotate_figure(F7A, fig.lab.size = 18, fig.lab = 'a'),
          annotate_figure(F7B+theme(legend.position = 'top'), fig.lab.size = 18, fig.lab = 'b'), widths = c(0.35,0.65)),
ggarrange(
  annotate_figure(FS13B, fig.lab.size = 18, fig.lab = 'c'),
  annotate_figure(FS13E, fig.lab.size = 18, fig.lab = 'd'),
  annotate_figure(FS13H, fig.lab.size = 18, fig.lab = 'e'),
  annotate_figure(F7F, fig.lab.size = 18, fig.lab = 'f'),
  annotate_figure(F7G, fig.lab.size = 18, fig.lab = 'g'),
  annotate_figure(F7H, fig.lab.size = 18, fig.lab = 'h'),
          nrow = 2, ncol = 3, heights = c(2,1)),
nrow = 2, heights = c(1,2.75))

#FS13
ggarrange(annotate_figure(FS13A,fig.lab.size = 18, fig.lab = 'a'),
          annotate_figure(FS13B,fig.lab.size = 18, fig.lab = 'b'),
          annotate_figure(FS13C,fig.lab.size = 18, fig.lab = 'c'),
          annotate_figure(FS13D,fig.lab.size = 18, fig.lab = 'd'),
          annotate_figure(FS13E,fig.lab.size = 18, fig.lab = 'e'),
          annotate_figure(FS13F,fig.lab.size = 18, fig.lab = 'f'),
          annotate_figure(FS13G,fig.lab.size = 18, fig.lab = 'g'),
          annotate_figure(FS13H,fig.lab.size = 18, fig.lab = 'h'),
          annotate_figure(FS13I,fig.lab.size = 18, fig.lab = 'i'),
          annotate_figure(FS13J,fig.lab.size = 18, fig.lab = 'j'),
          annotate_figure(FS13K,fig.lab.size = 18, fig.lab = 'k'),
          annotate_figure(FS13L,fig.lab.size = 18, fig.lab = 'l'),
nrow = 4, ncol = 3)
