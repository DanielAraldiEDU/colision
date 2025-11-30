# Colision

## Install FreeGLUT

```bash
brew install freeglut
```

## How to execute?

### Option 1: Using the VSCode

- Create the `tasks.json` file in `.vscode/`. Add this configuration it.

```json
{
  "version": "2.0.0",
  "tasks": [
    {
      "label": "Compilar com FreeGLUT",
      "type": "shell",
      "command": "g++",
      "args": [
        "main.cpp",
        "-o",
        "main",
        "-framework",
        "OpenGL",
        "-framework",
        "GLUT"
      ],
      "group": {
        "kind": "build",
        "isDefault": true
      }
    }
  ]
}
```

- Press `CMD + SHIFT + B` to compile.
- Run this command in terminal:

```bash
./main
```

### Option 2: Using the Makefile

Create the `Makefile` file with following command:

```makefile
CXX = g++
CXXFLAGS = -Wall -g
LDFLAGS = -framework OpenGL -framework GLUT

all: main

main: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o main $(LDFLAGS)

clean:
	rm -f main
```

- Compile with:

```bash
make
```

- Execute the program:

```bash
./main
```

- For clean program:

```bash
make clean
```
