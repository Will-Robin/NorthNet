ROOT=../NorthNet

CLASSES=$(ROOT)/Classes/*.py
LOADING=$(ROOT)/Loading/*.py $(ROOT)/Loading/network/*.py \
		$(ROOT)/Loading/reaction/*.py $(ROOT)/Loading/compound/*.py \
		$(ROOT)/Loading/experimentconditions/*.py
UTILS=$(ROOT)/Utils/*.py
WRITING=$(ROOT)/Writing/*.py $(ROOT)/Writing/model/*.py
MOL_OPS=$(ROOT)/molecule_operations/*.py
NET_GEN=$(ROOT)/network_generation/*.py
NETVIS=$(ROOT)/network_visualisation/*.py
RXN_OPS=$(ROOT)/reaction_operations/*.py
SIM=$(ROOT)/simulation/*.py

NORTHNET=$(CLASSES) $(LOADING) $(UTILS) $(WRITING) $(MOL_OPS) $(NET_GEN) \
		 $(NETVIS) $(RXN_OPS) $(SIM) $(ROOT)

all: index.html

index.html: $(NORTHNET)
	pdoc ../NorthNet -d "numpy" -o ./
