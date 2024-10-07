rm(list=ls())

require(MASS)
require(MBA)
require(gstat)
require(sf)
require(geoR)
require(ggplot2)

data("wolfcamp")
#class(wolfcamp)
#summary(wolfcamp)
#plot(wolfcamp)



wolfcamp_sf <- st_as_sf(data.frame(X1 = wolfcamp$coords[, 1], X2 = wolfcamp$coords[, 2],
                                   values = wolfcamp$data), coords = c("X1", "X2"))


idw_loocv <- function(data, p) {
  n <- nrow(data)  
  errors <- numeric(n)
  
  for (i in 1:n) {
    # Split the data: training set (without the i-th observation)
    train_data <- data[-i, ]
    test_data <- data[i, , drop = FALSE]  # drop = False, so that it does not collapse into a vector/point
    
    idw_model <- gstat(
      formula = values ~ 1, #~1 because we do not have any auxiliary variables
      locations = train_data, 
      set = list(idp = p)  #sets p
    )
    

    idw_prediction <- predict(idw_model, newdata = test_data)
    errors[i] <- ((idw_prediction$var1.pred - test_data$values) / sd(train_data$values))^2
  }
  
  return(mean(errors))
}


p_values <- seq(0.2, 10, by = 0.2)


cv_scores <- sapply(p_values, function(p) idw_loocv(wolfcamp_sf, p))


optimal_p <- p_values[which.min(cv_scores)]
optimal_cv <- min(cv_scores)



cv_data <- data.frame(
  p_values = p_values,  
  cv_scores = cv_scores 
)

ggplot(cv_data, aes(x = p_values, y = cv_scores)) +             
  geom_point(color = "blue", size = 2) +           
  geom_point(aes(x = optimal_p, y = optimal_cv), color = "red", size = 3) +
  labs(title = "Cross-Validation Score for Different p Values",
       x = "p",
       y = "CV") +
  theme_minimal()  

