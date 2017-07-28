orig_cols <- hsv(c(2,1,3)/3,
                 0.7 + 0.3 * (0:2/2),
                 0.3 + 0.7 * (2:0/2))
cols <- rep(orig_cols, each = 2)
ltys <- rep(c(2, 3), 3)
nms <- c("lin2", "lin2q", "lin3", "lin3q", "lin4", "lin4q")