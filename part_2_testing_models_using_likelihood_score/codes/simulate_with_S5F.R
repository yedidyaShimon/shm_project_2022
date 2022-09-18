library(shazam)
library(seqinr)

train_data_3<-read.csv("train_data_3_v_gene_final.csv",header=TRUE)

test_data_3_with_mut_num<-read.csv("test_data_3_v_gene_final.csv",header=TRUE)

#filtering
train_data_3<-train_data_3[train_data_3$num_mutations<20,]
test_data_3_with_mut_num<-test_data_3_with_mut_num[test_data_3_with_mut_num$num_mutations<20,]
print("filter")

tar_model_s<-createTargetingModel(train_data_3,model = "s",germlineColumn="ancestor_alignment",sequenceColumn ="sequence_alignment")

save(tar_model_s, file = "simulation_v_gene_3_object_s.RData")

print("model")

simData<-test_data_3_with_mut_num
simData$simulate_alignment<-rep("",nrow(simData))


for(i in 1:nrow(simData))
{
  simData$simulate_alignment[i]<-shmulateSeq(simData$ancestor_alignment[i],numMutations=simData$num_mutations[i],targetingModel=tar_model_s)
}


write.csv(simData, "test_data_3_S5F_s_simulate_final.csv", row.names=FALSE)
print("sim end")
