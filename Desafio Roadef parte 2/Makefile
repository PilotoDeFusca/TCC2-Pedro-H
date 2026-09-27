CXX        := g++
CXXFLAGS   := -std=c++17 -O2 -Wall -Wextra -I.
LDFLAGS    :=

TARGET     := solver
SRCS       := main.cpp loader.cpp trivial_solution.cpp verifier.cpp
DEPS       := loader.h

BUILD_DIR  := build
OBJS       := $(addprefix $(BUILD_DIR)/, $(SRCS:.cpp=.o))

all: $(BUILD_DIR)/$(TARGET)

$(BUILD_DIR)/$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: %.cpp $(DEPS)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: all
	./$(BUILD_DIR)/$(TARGET)

clean:
	rm -rf $(BUILD_DIR)

rebuild: clean all

.PHONY: all run clean rebuild