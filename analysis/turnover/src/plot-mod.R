
plot.all <- function(){
  ## f <- function(){
  ##   plot(mod, graph)
  ## }
  ## pdf.f(f, file= file.path(path, sprintf("%s.pdf",
  ##            paste("mod", "plot", sep="_"))),
  ##       width=20, height=20)

  f <- function(){
    dendPlot(mod, cex=0.5)
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste("dend", "plot", sep="_"))),
        width=20, height=20)

  f <- function(){
    plot(lc, type = "graph", layout = layout.fruchterman.reingold)
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste("mod1", "plot", sep="_"))),
        width=20, height=20)

  f <- function(){
    plot(lc, type = "graph", layout = "spencer.circle")
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste("mod2", "plot", sep="_"))),
        width=20, height=20)

  f <- function(){
    plot(lc, type = "graph", layout = "spencer.circle", shownodesin = 3)
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste("mod3", "plot", sep="_"))),
        width=20, height=20)

  f <- function(){
    plot(lc, type = "graph", shownodesin = 2, node.pies = TRUE)
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste("mod4", "plot", sep="_"))),
        width=20, height=20)

  f <- function(){
    plot(lc, type = "members")
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste("mod5", "plot", sep="_"))),
        width=20, height=20)

  f <- function(){
    plot(lc, type = "dend")
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste("mod6", "plot", sep="_"))),
        width=20, height=20)

  f <- function(){
    cr <- getClusterRelatedness(lc, hcmethod = "ward")
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste("dend2", "plot", sep="_"))),
        width=20, height=20)

  f <- function(){
    plot(lc, type = "commsumm", summary = "modularity")
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste("mod7", "plot", sep="_"))),
        width=20, height=20)

}

