setwd("~/data/2")
Sub_Merged_Mut_s<-read.table("Sub_Merged_Mut_s_all_t.tsv",header = TRUE)
Sub_Merged_Mut_rs<-read.table("Sub_Merged_Mut_rs_all_t.tsv",header = TRUE)

#setting the parameteres for the filtering
minNum_S<-200
minRatio_RS_to_S<-6


print(paste("requiring ",(minNum_S)," mutations, and ",minRatio_RS_to_S," times larger amount of silent + non-silent mutations than silent alone. Number of 5-mers that comply with requirements: ",sep=""))


Sub_Merged_Mut_s_measured<-Sub_Merged_Mut_s[Sub_Merged_Mut_s$fivemer.every & Sub_Merged_Mut_s$fivemer.total>minNum_S & Sub_Merged_Mut_rs$fivemer.total>minRatio_RS_to_S*Sub_Merged_Mut_s$fivemer.total,]
Sub_Merged_Mut_rs_measured<-Sub_Merged_Mut_rs[Sub_Merged_Mut_s$fivemer.every & Sub_Merged_Mut_s$fivemer.total>minNum_S & Sub_Merged_Mut_rs$fivemer.total>minRatio_RS_to_S*Sub_Merged_Mut_s$fivemer.total,]
print(nrow(Sub_Merged_Mut_s_measured))
print(nrow(Sub_Merged_Mut_rs_measured))
list_s_A<-Sub_Merged_Mut_s_measured[,1:4]$A
list_s_A<-list_s_A[complete.cases(list_s_A)]

list_rs_A<-Sub_Merged_Mut_rs_measured[,1:4]$A
list_rs_A<-list_rs_A[complete.cases(list_rs_A)]

list_s_T<-Sub_Merged_Mut_s_measured[,1:4]$T
list_s_T<-list_s_T[complete.cases(list_s_T)]

list_rs_T<-Sub_Merged_Mut_rs_measured[,1:4]$T
list_rs_T<-list_rs_T[complete.cases(list_rs_T)]

list_s_C<-Sub_Merged_Mut_s_measured[,1:4]$C
list_s_C<-list_s_C[complete.cases(list_s_C)]

list_rs_C<-Sub_Merged_Mut_rs_measured[,1:4]$C
list_rs_C<-list_rs_C[complete.cases(list_rs_C)]

list_s_G<-Sub_Merged_Mut_s_measured[,1:4]$G
list_s_G<-list_s_G[complete.cases(list_s_G)]

list_rs_G<-Sub_Merged_Mut_rs_measured[,1:4]$G
list_rs_G<-list_rs_G[complete.cases(list_rs_G)]

list_s<-c(list_rs_A,list_s_C,list_s_G,list_s_T)
list_rs<-c(list_rs_A,list_rs_C,list_rs_G,list_rs_T)
print("The correlation coefficient of the substitution to A rates is:")
print(cor(list_rs_A,list_s_A))
print("The correlation coefficient of the substitution to C rates is:")
print(cor(list_rs_C,list_s_C))
print("The correlation coefficient of the substitution to G rates is:")
print(cor(list_rs_G,list_s_G))
print("The correlation coefficient of the substitution to T rates is:")
print(cor(list_rs_T,list_s_T))
print("The general correlation coefficient rate is")
print(cor(list_rs,list_s))






#cor_s_and_sr<-cancor((,(Sub_Merged_Mut_rs_measured[,1:4])$A)

#c(Sub_Merged_Mut_s_measured[,1:4])
