# Development principles

## Architecture

Spoc aims to separate domain logic from adaptors (like file system or database interactions) based on the principles of [hexaongal software architecutre](https://en.wikipedia.org/wiki/Hexagonal_architecture_(software)). For file-system interactions, this is implemented by having a single point of IO in the [io.pyfile](../spoc/io.py) that implements reading, writing and validating files from disk. The application layer of spoc that deals with domain-logic is separated from file-system interactions and only receives and returns python objects. Whenever possible, the domain logic classes should follow functional principles and try not perform their tasks using side-effects. This facilitates testability of spoc.