library(venneuler)

A <- c(1:10)
B <- c(4:14)
C <- c(13:15)

sum(A %in% B)
sum(B %in% A)

venn <- venneuler(c(A = length(A),
                    B = length(B),
                    C = length(C),
                    "A&B" = sum(A%in%B),
                    "A&C" = sum(A%in%C),
                    "B&C" = sum(B%in%C)))

plot(venn)
