
# Name of model/process
NAME = Root
MODELPATH = "Root"
SOLIB = usrLib$(NAME).so

# List all object files here
OBJECTS = $(NAME).o PBD.o tissue.o
# default target
target: $(SOLIB)

# Call MorphoDynamX to get makefile
include $(shell mdx --resource)/MDXProcess.mk


# Add extra compile flags here
CXXFLAGS += -Wno-deprecated-copy -Wno-unused-local-typedefs -Wno-unused-parameter $(EXTRA) $(INC)


# Add extra libraries here
#LIBS+=Root.cu.o -lcudart

# Add extra link flags here
LD_FLAGS+=   $(EXTRA)

# Model dependencies
$(NAME).o: $(NAME).cpp $(NAME).hpp  Makefile

$(SOLIB): $(OBJECTS)

# Run the model
run: target
	mdx --model 'Model/Root/01 Root' --addlibrary $(SOLIB) RootModel.mdxv
