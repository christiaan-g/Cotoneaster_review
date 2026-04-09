nn_dist_to_arc<-
  function (x, r = 6371000){
    d <- acos( 1 - ( (x^2) / ( 2*(r^2)) ) ) *r
    round(d)
  }

# # Validating the formula
# 
# nn_dist_to_arc(10)
# nn_dist_to_arc(x = 10000 *10^3)
# 
# result <- data.frame(
#   id = seq(2,2*63711,1),
#   dist = seq(0,2*6371000,100),
#   arc_dist = NA
# )
# 
# result$arc_dist <- nn_dist_to_arc(result$dist)
# 
# plot(result$arc_dist,result$dist)
# 
# result$arc_dist[nrow(result)] / 6371000  # should be close to pi
