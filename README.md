# NQueens-Genetic-Algorithm

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
- [Usage](#usage)
- [Components](#components)
- [Customization](#customization)
- [License](#license)

## Introduction

The **NQueens-Genetic-Algorithm** is an intelligent and efficient solution to the classic N-Queens problem. The N-Queens problem involves placing N chess queens on an N×N chessboard, ensuring that no two queens threaten each other. Solving this problem is fundamental to artificial intelligence and has real-world applications in constraint satisfaction problems.

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
   git clone https://github.com/Michael-Moulana/NQueens-Genetic-Algorithm.git

   ```

2. Open the project in your development environment and configure the parameters in the genetic algorithm components to meet your requirements.

3. Compile and run the program to solve the N-Queens problem for your chosen N value.

## Usage

### Running the Program:

- Execute the program to initiate the genetic algorithm.
- The algorithm will evolve a population of candidate solutions over generations, aiming to find a valid N-Queens solution with no threats between queens.

### Customization:

- You can fine-tune the genetic algorithm by modifying parameters in the components to better suit your specific problem size and requirements.

## Components

The **NQueens-Genetic-Algorithm** comprises the following core components:

- **RouletteWheelSelection**: Implements roulette wheel selection for choosing chromosomes based on their fitness, with built-in selective pressure control.

- **TournamentSelection**: Utilizes tournament selection for chromosome selection, enabling diverse and robust solutions.

- **Repair**: Ensures that chromosomes are free from duplicate values by using a sophisticated repair mechanism that adheres to the problem's constraints.

- **Mutation**: Implements mutation operators that introduce subtle changes to the population to encourage exploration and avoid premature convergence.

- **Fitness**: Evaluates the fitness of a chromosome by counting the number of attacks in the N-Queens problem, a critical measure of solution quality.

## Customization

One of the key strengths of this program is its customizability. You can adjust various parameters such as population size, mutation rate, selective pressure, and crossover methods. Experimenting with these parameters can lead to improved results and tailored solutions for your N-Queens problem.

## License

This project is licensed under the MIT License. See the LICENSE.md file for more details.
