#script for building volcano-plot using ggplot2
#Author: Benjamin Clark

#License

# © Copyright 2020 iGEM Concordia, Maher Hassanain, Benjamin Clark, Hajar Mouddene, Grecia Orlano
# This file is part of AstroBio.
# AstroBio is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
# AstroBio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with AstroBio.  If not, see <https://www.gnu.org/licenses/>.

library(ggplot2)
library(dplyr)



build_volcano <- function(rem_data){
  data <- rem_data %>% select(randomSummary, randomP, signcon, Symbol) 
  data$randomP <- -log10(data$randomP)
  
  p <- ggplot(data = data, aes(x = randomSummary, y = randomP, col = signcon)) + geom_point() + theme_minimal() + 
    scale_color_gradient2(mid = "grey95", low = "darkblue", high = "gold3") +
    labs(colour =  "Sign Consistency", x = "Summary Log2 Fold-Change", y = "-Log10 Summary P-Value")
  
  return(p)
} 



