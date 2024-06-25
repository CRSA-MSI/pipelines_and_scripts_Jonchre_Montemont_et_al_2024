library("stats")

#### #### #### ######## #### MSI/N ######## #### #### ######## #### #### ####
### EXclusion 
#### 3'MICE vs 5'Mice
chisq.test(matrix(c(513,79661,44,54739),2,2, byrow=TRUE), correct=TRUE) 
#p-value < 2.2e-16

### EXclusion 
#### 5'Mice vs 5'3'MICE
chisq.test(matrix(c(44,54739, 326, 47355),2,2, byrow=TRUE), correct=TRUE) 
#p-value < 2.2e-16

### EXclusion 
#### 3'MICE vs 5'3'MICE
chisq.test(matrix(c(513,79661,326, 47355),2,2, byrow=TRUE), correct=TRUE) 
#p-value = 0.3664

#### #### #### ####
### Inclusion 
#### 3'MICE vs 5'Mice
chisq.test(matrix(c(84,79661,6,54739),2,2, byrow=TRUE), correct=TRUE) 
#p-value = 9.938e-11

### Inclusion 
#### 5'Mice vs 5'3'MICE
chisq.test(matrix(c(6,54739, 58, 47355),2,2, byrow=TRUE), correct=TRUE) 
#p-value = 3.185e-12

### Inclusion 
#### 3'MICE vs 5'3'MICE
chisq.test(matrix(c(84,79661,58, 47355),2,2, byrow=TRUE), correct=TRUE) 
#p-value = 0.4292

#### #### #### ######## #### MSS/N ######## #### #### ######## #### #### ####
### EXclusion 
#### 3'MICE vs 5'Mice
chisq.test(matrix(c(75,79661,13,54739),2,2, byrow=TRUE), correct=TRUE) 
#p-value = 1.261e-06

### EXclusion 
#### 5'Mice vs 5'3'MICE
chisq.test(matrix(c(13,54739, 39, 47355),2,2, byrow=TRUE), correct=TRUE) 
#  p-value = 6.395e-05
 
### EXclusion 
#### 3'MICE vs 5'3'MICE
chisq.test(matrix(c(75,79661,39, 47355),2,2, byrow=TRUE), correct=TRUE) 
#p-value = 0.5611

#### #### #### ####
### Inclusion 
#### 3'MICE vs 5'Mice
chisq.test(matrix(c(8,79661,2,54739),2,2, byrow=TRUE), correct=TRUE) 
#p-value = 0.3114

### Inclusion 
#### 5'Mice vs 5'3'MICE
chisq.test(matrix(c(2,54739, 8, 47355),2,2, byrow=TRUE), correct=TRUE) 
#p-value = 0.0696

### Inclusion 
#### 3'MICE vs 5'3'MICE
chisq.test(matrix(c(8,79661,8, 47355),2,2, byrow=TRUE), correct=TRUE) 
#p-value = 0.4276

#### #### #### ######## #### MSI/MSS ######## #### #### ######## #### #### ####
### EXclusion 
#### 3'MICE vs 5'Mice
chisq.test(matrix(c(154,79661,7,54739),2,2, byrow=TRUE), correct=TRUE) 
# p-value < 2.2e-16

### EXclusion 
#### 5'Mice vs 5'3'MICE
chisq.test(matrix(c(7,54739, 122, 47355),2,2, byrow=TRUE), correct=TRUE) 
# p-value < 2.2e-16

### EXclusion 
#### 3'MICE vs 5'3'MICE
chisq.test(matrix(c(154,79661,122, 47355),2,2, byrow=TRUE), correct=TRUE) 
#p-value = 0.02075
 
#### #### #### ####
### Inclusion 
#### 3'MICE vs 5'Mice
chisq.test(matrix(c(8,79661,2,54739),2,2, byrow=TRUE), correct=TRUE) 
#p-value = 0.3114
 
### Inclusion 
#### 5'Mice vs 5'3'MICE
chisq.test(matrix(c(2,54739, 12, 47355),2,2, byrow=TRUE), correct=TRUE) 
#p-value = 0.007302
 
### Inclusion 
#### 3'MICE vs 5'3'MICE
chisq.test(matrix(c(8,79661,12, 47355),2,2, byrow=TRUE), correct=TRUE) 
# p-value = 0.06154
 
MSI <- c(6.44,1.05,0.80,0.11,6.88,1.22)/100
names(MSI) <- c("3.MICE-Exclusion-MSI/N","3.MICE-Inclusion-MSI/N","5.Mice-Exclusion-MSI/N",
		"5.Mice-Inclusion-MSI/N","5.-3.MICE-Exclusion-MSI/N","5.-3.MICE-Inclusion-MSI/N")					
MSS <- c(0.94,0.10,0.24,0.04,0.82,0.17)/100
names(MSS) <- c("3.MICE-Exclusion-MSS/N","3.MICE-Inclusion-MSS/N","5.Mice-Exclusion-MSS/N",
		"5.Mice-Inclusion-MSS/N","5.-3. MICE-Exclusion-MSS/N","5.-3. MICE-Inclusion-MSS/N") 

svg("Fig.svg")
par(mfrow = c(2,2))
par(oma=c(5,1,1,0))
barplot(MSI[c(1,3,5)], las=2, ylab=("Percentage of exons"), main="MSI.EXCLUSION")
barplot(MSI[c(2,4,6)], las=2, ylab=("Percentage of exons"), main="MSI.INCLUSION")
barplot(MSS[c(1,3,5)], las=2, ylab=("Percentage of exons"), main="MSS.EXCLUSION")
barplot(MSS[c(2,4,6)], las=2, ylab=("Percentage of exons"), main="MSS.INCLUSION")
dev.off()

