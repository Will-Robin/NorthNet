ROOT=../NorthNet

CLASSES=$(ROOT)/Classes/*.py
#LOADING=$(ROOT)/Loading/*.py
UTILS=$(ROOT)/Utils/*.py
WRITING=$(ROOT)/Writing/*.py $(ROOT)/Writing/model/*.py
MOL_OPS=$(ROOT)/molecule_operations/*.py
NET_GEN=$(ROOT)/network_generation/*.py
NETVIS=$(ROOT)/network_visualisation/*.py
RXN_OPS=$(ROOT)/reaction_operations/*.py
TXT_PARSE=$(ROOT)/text_parsing/*.py
SIM=$(ROOT)/simulation/*.py

NORTHNET=$(CLASSES) $(LOADING) $(UTILS) $(WRITING) $(MOL_OPS) $(NET_GEN) \
		 $(NETVIS) $(RXN_OPS) $(TXT_PARSE) $(SIM)

all: index.html

index.html: $(NORTHNET)
	pdoc ../NorthNet -d "numpy" -o ./
