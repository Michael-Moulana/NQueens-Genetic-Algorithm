# Genetic Algorithm for N-Queens Problem Solver

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
- [Usage](#usage)
- [Components](#components)
- [Customization](#customization)
- [Contributing](#contributing)
- [License](#license)

## Introduction

The **Genetic Algorithm for N-Queens Problem Solver** is an intelligent and efficient solution to the classic N-Queens problem. The N-Queens problem involves placing N chess queens on an N×N chessboard, ensuring that no two queens threaten each other. Solving this problem is fundamental to artificial intelligence and has real-world applications in constraint satisfaction problems.

Our program employs a genetic algorithm approach to tackle the N-Queens problem, offering an elegant and scalable solution. Genetic algorithms are inspired by the process of natural selection and evolution, making them highly adaptable to complex optimization challenges like this one.

## Features

- **Efficient**: This solution is optimized for solving N-Queens problems of varying sizes (N) with speed and accuracy.
- **Crossover Methods**: Supports a variety of crossover methods, including single-point, uniform, and order1, allowing you to experiment with different strategies.
- **Mutation Operators**: Apply mutation to introduce small variations in the population, driving exploration and diversification.
- **Repair Mechanism**: Ensures the integrity of the chromosome solutions by eliminating duplicate values and preserving the problem's constraints.
- **Customization**: Tailor the genetic algorithm's parameters, population size, mutation rate, and selection strategies to fit your specific problem or experiment.
- **Easy Integration**: The well-organized C# codebase and clear class separation make it straightforward to integrate into your AI projects.

## Getting Started

### Prerequisites

- [Visual Studio](https://visualstudio.microsoft.com/) or your preferred C# development environment.
- Basic knowledge of genetic algorithms and the N-Queens problem.

### Installation

1. Clone the repository to your local machine:

   ```shell
   git clone https://github.com/yourusername/genetic-nqueens.git

2.Open the project in your development environment and configure the parameters in the genetic algorithm components to meet your requirements.

3.Compile and run the program to solve the N-Queens problem for your chosen N value.

## Usage
Running the Program:

Execute the program to initiate the genetic algorithm.
The algorithm will evolve a population of candidate solutions over generations, aiming to find a valid N-Queens solution with no threats between queens.
Customization:

You can fine-tune the genetic algorithm by modifying parameters in the components to better suit your specific problem size and requirements.
