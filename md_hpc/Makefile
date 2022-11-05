STRUTS=../generate_mutation/Structures
MUTREC=../generate_mutation/Recipes
IMAGE_TAG=gab0/gromacs-toxo-hpc:v28
REFERENCE_STRUCT=/home/Database/Research/TGondiiAntigens/Protocolo/generate_mutation/Structures/original_0_SAG1__0001.pdb
F1=T_400_SALT.obj
F2=T_4000_SALT.obj
all:
	ls
image:
	docker build -t gromacs-hpc .
upload-image:
	docker tag gromacs-hpc $(IMAGE_TAG)
	docker push $(IMAGE_TAG)
sif:
	SINGULARITY_NOHTTPS=1 singularity build --sandbox gromacs-locald.sif docker-daemon://gromacs-local:latest
	sudo singularity build gromacs-local.sif gromacs-locald.sif
retrieve:
	rsync --rsh='ssh -p31459' --append-verify -av gabzn@cenapad.unicamp.br:~/homelovelace/MD/results .

define analyze
	mdanalyze -i $(3) \
			--load-session results/$(1) \
			--reference-structure $(REFERENCE_STRUCT) \
			--selection $(2)
endef

generate:
# This part is project-specific,
# matching the original project's structure groups.
	$(call analyze,$(F1),1-4,DERIV)
	$(call analyze,$(F1),17-20,NAT)
	$(call analyze,$(F1),"8;5-7",D)
	$(call analyze,$(F1),9-16,DALL)
	$(call analyze,$(F2),all,NAT4K)
