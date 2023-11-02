using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace GeneticAlgorithm
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.SetWindowSize(120, 13);
            Console.ForegroundColor = ConsoleColor.White;
            Console.BackgroundColor = ConsoleColor.DarkGray;
            Console.Clear();

            int experimentsQuantity = 10; // Number of experiments - default: 10

            int numberOfQueens = 50; // Number of queens - default: 100
            int generationPopulation = 500; // number of members in each generation - default: 500
            int numberOfGenerations = 100; // number of generations that we perform genetic algorithm on it - default: 100

            double alphaValue = 1.2; // specifies the alpha value in roulette-wheel selection - default: 1.4
            double betaValue = 15; // specifies the beta value in roulette-wheel selection - default: 10
            int k = 400; // specifies the "k" value in k-tournament selection - default: 10

            double crossoverRate = 0.6; // Rate of Crossover - default: 0.5
            double mutationRate = 0.7; // Rate of Mutation - default: 0.1

            //string encoding = "standard"; // 2 modes: "standard" & "permutation"
            string encoding = "permutation"; // 2 modes: "standard" & "permutation"

            string selectionMethod = "tournament"; // 2 methods: "rouletteWheel" & "tournament"
            //string selectionMethod = "rouletteWheel"; // 2 methods: "rouletteWheel" & "tournament"

            //string crossoverMethod = "singlePoint"; // 3 methods: "singlePoint" & "uniform" & "order1"
            //string crossoverMethod = "uniform"; // 3 methods: "singlePoint" & "uniform" & "order1"
            string crossoverMethod = "order1"; // 3 methods: "singlePoint" & "uniform" & "order1"

            List<GeneticAlgorithm> GA_Experiments = new List<GeneticAlgorithm>();

            DateTime startingTime1 = DateTime.UtcNow; // starting the timer
            for (int i = 0; i < experimentsQuantity; i++)
            {
                Console.ForegroundColor = ConsoleColor.White;
                Console.WriteLine("                                              Running Experiment: {0}/{1}", i + 1, experimentsQuantity);

                GA_Experiments.Add(new GeneticAlgorithm(
                 numberOfQueens, generationPopulation, numberOfGenerations,
                 alphaValue, betaValue, k,
                 crossoverRate, mutationRate,
                 encoding, selectionMethod, crossoverMethod));

                Console.ForegroundColor = ConsoleColor.Black;
                GA_Experiments[i].PerformAlgorithm();

                Console.Clear();
            }
            DateTime endTime1 = DateTime.UtcNow; // stopping the timer

            TimeSpan timeDifference1 = endTime1 - startingTime1;
            double seconds = timeDifference1.TotalSeconds; // total elapsed time in seconds 

            Result.Analyzer(GA_Experiments, experimentsQuantity, seconds);

            Console.ReadKey();
        }
    }

    class Result
    {
        static int experimentsQuantity;

        static int numberOfQueens; // Number of queens - default: 100
        static int generationPopulation; // number of members in each generation - default: 500
        static int numberOfGenerations; // number of generations that we perform genetic algorithm on it - default: 100

        static string alphaValue; // specifies the alpha value in roulette-wheel selection
        static string betaValue ; // specifies the beta value in roulette-wheel selection
        static string k; // specifies the "k" value in k-tournament selection

        static double crossoverRate; // Rate of Crossover - default: 0.5
        static double mutationRate; // Rate of Mutation - default: 0.1

        static string encoding; // 2 modes: "standard" & "permutation"
        static string selectionMethod; // 2 methods: "rouletteWheel" & "tournament"
        static string crossoverMethod; // 3 methods: "singlePoint" & "uniform" & "order1"
        static string mutationMethod; // 2 methods: "standard" & "swap"

        static List<double> avgInitialFitness = new List<double>(); // average of "initial fitness values" for all generations
        static List<double> avgFitness = new List<double>(); // average of "fitness values" for all generations
        static List<double> avgMinFitness = new List<double>(); // average of "last fitness values" for all generations
        static List<int> successValues = new List<int>(); // determines the success rate

        public static void Analyzer(List<GeneticAlgorithm> results, int numberOfExperiments, double elapsedTime)
        {
            experimentsQuantity = numberOfExperiments;

            numberOfQueens = results[0].numberOfQueens;
            generationPopulation = results[0].generationPopulation;
            numberOfGenerations = results[0].numberOfGenerations;

            alphaValue = Convert.ToString(results[0].alphaValue);
            betaValue = Convert.ToString(results[0].betaValue);
            k = Convert.ToString(results[0].k);

            crossoverRate = results[0].crossoverRate;
            mutationRate = results[0].mutationRate;

            encoding = results[0].encoding;

            selectionMethod = results[0].selectionMethod;
            if (selectionMethod == "rouletteWheel")
                k = "Not Configured";
            else if(selectionMethod == "tournament")
            {
                alphaValue = "Not Configured";
                betaValue = "Not Configured";
            }

            crossoverMethod = results[0].crossoverMethod;
            if (encoding == "standard")
                mutationMethod = "standard";
            else if(encoding == "permutation");
                mutationMethod = "swap";

            for (int i = 0; i < numberOfExperiments; i++)
            {
                avgInitialFitness.Add(results[i].initialFitness);
                avgFitness.Add(results[i].avgFitness);
                avgMinFitness.Add(results[i].lastFitness);
                successValues.Add(results[i].done);
            }

            Console.ForegroundColor = ConsoleColor.White;

            Console.WriteLine(" ******************    Genetic Algorithm Performance Summary: (Performed on {0}-Queens Problem)    *********************",
                numberOfQueens);

            Console.ForegroundColor = ConsoleColor.Black;

            Console.WriteLine("\n                 Number of Experiments: {0} | Generation Population: {1} | Number of Generations: {2}  ",
                experimentsQuantity, generationPopulation, numberOfGenerations);
            Console.WriteLine("\n                 Encoding: {0} | Selection: {1} | Crossover: {2} | mutation: {3}     ",
                encoding, selectionMethod, crossoverMethod, mutationMethod);
            Console.WriteLine("\n     Alpha Value: {0} | Beta Value: {1} | K: {2} | Crossover Rate: {3}% | Mutation Rate: {4}%      ",
                alphaValue, betaValue, k, crossoverRate * 100, mutationRate * 100);
            Console.WriteLine("\n           Success Rate: {0}% | Execution Time for each GA: {1} Seconds | Total Execution Time: {2} Seconds       ",
                successValues.Average() * 100, (elapsedTime/experimentsQuantity).ToString("N2"), elapsedTime.ToString("N1"));
            Console.WriteLine("\n                       Initial Fitness Avg: {0} | Fitness Avg: {1} | Best Fitness Avg: {2}       ",
                Math.Round(avgInitialFitness.Average()/generationPopulation, 2), Math.Round(avgFitness.Average()/numberOfGenerations, 2), avgMinFitness.Average());
        }
    }

    public class ProgressBar : IDisposable, IProgress<double>
    {
        private const int blockCount = 10;
        private readonly TimeSpan animationInterval = TimeSpan.FromSeconds(1.0 / 8);
        private const string animation = @"|/-\";

        private readonly Timer timer;

        private double currentProgress = 0;
        private string currentText = string.Empty;
        private bool disposed = false;
        private int animationIndex = 0;

        public ProgressBar()
        {
            timer = new Timer(TimerHandler);

            // A progress bar is only for temporary display in a console window.
            // If the console output is redirected to a file, draw nothing.
            // Otherwise, we'll end up with a lot of garbage in the target file.
            if (!Console.IsOutputRedirected)
            {
                ResetTimer();
            }
        }

        public void Report(double value)
        {
            // Make sure value is in [0..1] range
            value = Math.Max(0, Math.Min(1, value));
            Interlocked.Exchange(ref currentProgress, value);
        }

        private void TimerHandler(object state)
        {
            lock (timer)
            {
                if (disposed) return;

                int progressBlockCount = (int)(currentProgress * blockCount);
                int percent = (int)(currentProgress * 100);
                string text = string.Format("                                                 [{0}{1}] {2,3}% {3}",
                    new string('#', progressBlockCount), new string('-', blockCount - progressBlockCount),
                    percent,
                    animation[animationIndex++ % animation.Length]);
                UpdateText(text);

                ResetTimer();
            }
        }

        private void UpdateText(string text)
        {
            // Get length of common portion
            int commonPrefixLength = 0;
            int commonLength = Math.Min(currentText.Length, text.Length);
            while (commonPrefixLength < commonLength && text[commonPrefixLength] == currentText[commonPrefixLength])
            {
                commonPrefixLength++;
            }

            // Backtrack to the first differing character
            StringBuilder outputBuilder = new StringBuilder();
            outputBuilder.Append('\b', currentText.Length - commonPrefixLength);

            // Output new suffix
            outputBuilder.Append(text.Substring(commonPrefixLength));

            // If the new text is shorter than the old one: delete overlapping characters
            int overlapCount = currentText.Length - text.Length;
            if (overlapCount > 0)
            {
                outputBuilder.Append(' ', overlapCount);
                outputBuilder.Append('\b', overlapCount);
            }

            Console.Write(outputBuilder);
            currentText = text;
        }

        private void ResetTimer()
        {
            timer.Change(animationInterval, TimeSpan.FromMilliseconds(-1));
        }

        public void Dispose()
        {
            lock (timer)
            {
                disposed = true;
                UpdateText(string.Empty);
            }
        }

    }

    class GeneticAlgorithm
    {
        public int numberOfQueens; // Number of queens - default: 100
        public int generationPopulation; // number of members in each generation - default: 500
        public int numberOfGenerations; // number of generations that we perform genetic algorithm on them - default: 100

        public double alphaValue; // specifies the alpha value in roulette-wheel selection
        public double betaValue; // specifies the beta value in roulette-wheel selection
        public int k; // specifies the "k" value in k-tournament selection

        public double crossoverRate; // Rate of Crossover - default: 0.5
        public double mutationRate; // Rate of Mutation - default: 0.1

        public Generation currentGeneration; // list of chromosomes together, which is called a generation

        public string encoding; // 2 modes: "standard" & "permutation"
        public string selectionMethod; // 2 methods: "rouletteWheel" & "tournament"
        public string crossoverMethod; // 3 methods: "singlePoint" & "uniform" & "order1"

        public static Chromosome[] modifiedChromosomes; // storing chromosomes that crossover applied on them, temporarily

        // public int fitnessFuncCounter = 0; // number of fitness function executions
        public double initialFitness; // the initial fitness value for a given generation
        // public double lastFitness; // the last "best" fitness value in the last generation
        public int lastFitness = int.MaxValue; // the best fitness result in a chromosome whitin a generation

        public int totalFitness = 0; // sum of all fitness values in a generation
        public int avgFitness = 0; // avg of all fitness values in a generation 

        public int done = 0; // in case of algorithm completion, this value changes to "true", otherwise "false"

        public double probability; // probability of doing crossover or mutation
        public static Random rnd = new Random();

        public GeneticAlgorithm(
            int numberOfQueens, int generationPopulation, int numberOfGenerations,
            double alphaValue, double betaValue, int k,
            double crossoverRate, double mutationRate,
            string encoding, string selectionMethod, string crossoverMethod) 
        {
            this.numberOfQueens = numberOfQueens;
            this.generationPopulation = generationPopulation;
            this.numberOfGenerations = numberOfGenerations;
            this.alphaValue = alphaValue;
            this.betaValue = betaValue;
            this.k = k;
            this.crossoverRate = crossoverRate;
            this.mutationRate = mutationRate;
            this.encoding = encoding;
            this.selectionMethod = selectionMethod;
            this.crossoverMethod = crossoverMethod;
            CreateGeneration();
        }

        public void CreateGeneration()
        {
            currentGeneration = new Generation(numberOfQueens, generationPopulation, encoding);
        }

        public void PerformAlgorithm()
        {
            using (var progress = new ProgressBar())
            {
                    for (int generationCounter = 0; generationCounter < numberOfGenerations; generationCounter++)
                    {

                    if (generationCounter % 5 == 0)
                    {
                        progress.Report((double)generationCounter / numberOfGenerations);
                        Thread.Sleep(20);
                    }

                    totalFitness = 0; // resetting the value in every generation
                        lastFitness = int.MaxValue;
                        for (int i = 0; i < generationPopulation; i++) // examining all of the chromosomes' fitness values in the last generation
                        {
                            if (generationCounter == 0)
                                initialFitness += currentGeneration.chromosomes[i].fitness;

                            //fitnessFuncCounter++;
                            totalFitness += currentGeneration.chromosomes[i].fitness; // sum of all fitness values in a generation

                            if (currentGeneration.chromosomes[i].fitness <= lastFitness) // finding the best fitness in a generation 
                                lastFitness = currentGeneration.chromosomes[i].fitness;

                            if (lastFitness == 0) // if an optimal solution found, terminate the loop
                            {
                                done = 1;
                                break;
                            }
                            //Console.WriteLine(i);
                        }
                        //Console.WriteLine("best fitness: {0}", lastFitness);
                        //Console.WriteLine(totalFitness);
                        //Console.WriteLine("/////////////////////////////////////////////////////////////////  {0}", generationCounter);
                        avgFitness = avgFitness + (totalFitness / generationPopulation); // avg fitness of all chromosomes in multiple generations

                        if (done == 1)
                            break;

                        currentGeneration = selection.ApplySelection(currentGeneration, alphaValue, betaValue, k, selectionMethod); // applying selection

                        for (int i = 0; i < generationPopulation; i = i + 2) // applying crossover
                        {
                            probability = rnd.NextDouble();
                            if (probability <= crossoverRate)
                            {
                                modifiedChromosomes = Crossover.ApplyCrossover(
                                    currentGeneration.chromosomes[i],
                                    currentGeneration.chromosomes[i + 1],
                                    crossoverMethod);

                                currentGeneration.chromosomes[i] = modifiedChromosomes[0];
                                currentGeneration.chromosomes[i + 1] = modifiedChromosomes[1];
                            }
                        }

                        for (int i = 0; i < generationPopulation; i++) // applying mutation
                        {
                            probability = rnd.NextDouble();
                            if (probability <= mutationRate)
                            {
                                currentGeneration.chromosomes[i] = Mutation.ApplyMutation(currentGeneration.chromosomes[i]);
                            }
                        }
                    }
            }


        }
    }

    class Chromosome
    {
        public int[] chromosome; // stores the state of queens in a one dimensional array - array's length will be determined by user inputs through the constructor

        // determines in which process this chromosome created
        // 1st process: "standard"
        // 2nd process: "permutation"
        public string creationProcess;

        public int fitness; // fitness of the chromosome, evaluated by fitness function


        public static Random rnd = new Random();

        // Constructing chromosomes with pre-determined queens' locations (i.e. after operators' modifications applied)
        public Chromosome(int[] newChromosome, string creationProcess)
        {
            chromosome = newChromosome;
            this.creationProcess = creationProcess;
            fitness = Fitness.FitnessFunction(newChromosome);
        }

        public Chromosome(int N, string creationProcess)
        {
            if (creationProcess == "standard")
                Standard(N);
            else
                Permutation(N);
        }

        public void Standard(int N) // This constructor will generate an standard randomized starting state, based on the "number of queens (N)"
        {
            chromosome = new int[N]; // Building The First Chromosome, Length Determined by N 

            // Building The First chromosome in a randomized order of numbers between "0 and N-1" "i.e {1,0,1,4,0,5,6,7} for N=8" "repetetive numbers may occur"
            for (int i = 0; i < N; i++)
                chromosome[i] = rnd.Next(N);

            creationProcess = "standard";
            fitness = Fitness.FitnessFunction(chromosome);
        }

        public void Permutation(int N) // This constructor will generate a randomized starting state with permutation, based on the "number of queens (N)"
        {
            chromosome = new int[N]; // Building The First Chromosome, Length Determined by N 
            for (int i = 0; i < N; i++) // Building The First chromosome in order of 0,...,N-1 "i.e {0,1,2,3,4,5,6,7} for N=8"
                chromosome[i] = i;

            chromosome = Enumerable.Range(0, N).OrderBy(r => rnd.Next()).ToArray(); // Randomizing The First State "i.e {2,3,5,7,4,1,6,0} for N=8"

            creationProcess = "permutation";
            fitness = Fitness.FitnessFunction(chromosome);
        }

        //public Chromosome(int N) // This constructor will generate a randomized starting state, based on the "number of queens (N)"
        //{
        //    // chromosome = new int[] { 6, 4, 2, 0, 5, 7, 1, 3 }; // 0 vazir guard
        //    chromosome = new int[] { 6, 0, 5, 2, 4, 3, 6, 1 }; // 6 barkhord
        //    // chromosome = new int[] { 0, 2, 4, 1, 3, 6, 7, 5 }; // 4 barkhord
        //}
    }

    class Generation
    {
        public int N; // Number of queens - default: 100
        public int generationPopulation; // number of members in each generation - default: 500
        public string encoding; // 2 modes: "standard" & "permutation"

        public List<Chromosome> chromosomes = new List<Chromosome>(); // list of chromosomes together, which is called a generation

        public Generation(int N, int generationPopulation, string encoding) // creating first generation
        {
            this.N = N;
            this.generationPopulation = generationPopulation;
            this.encoding = encoding;

            for (int i = 0; i < generationPopulation; i++)
            {
                chromosomes.Add(new Chromosome(N, encoding));
            }
        }

        public Generation(List<Chromosome> chromosomes) // assigning a new generation made by selection methods
        {
            N = chromosomes[0].chromosome.Length;
            generationPopulation = chromosomes.Count;
            encoding = chromosomes[0].creationProcess;

            this.chromosomes = chromosomes;
        }
    }

    class selection
    {
        public static Generation ApplySelection(Generation currentGeneration, double alphaValue, double betaValue, int k, string selectionMethod)
        {
            if (selectionMethod == "rouletteWheel")
            {
                return RouletteWheelSelection.ApplyRouletteWheelSelection(currentGeneration, alphaValue, betaValue);
            }
            else if(selectionMethod == "tournament")
            {
                return TournamentSelection.ApplyTournamentSelection(currentGeneration, k);
            }

            return null;
        }
    }

    class RouletteWheelSelection
    {
        public static int generationPopulation;
        // to apply selective pressure, default: 1.4 (determining "how much" better fitness values than the worst fitness in generation, gets prioritized)
        public static double alpha;
        // to apply selective pressure, default: 10 (if beta sets higher, algorithm will dedicate higher chance for better chromosomes)
        public static double beta;
        // the worst fitness function value in the current generation
        public static int currentGeneration_WorstFitness;

        public static double portion;
        public static double[] portionPercentage; // can be a fraction of range:(0...1) for each chromosome based on their fitness value

        public static List<Chromosome> winnerChromosomes;

        public static Random rnd = new Random();

        public static Generation ApplyRouletteWheelSelection(Generation currentGeneration, double alphaValue, double betaValue)
        {
            generationPopulation = currentGeneration.chromosomes.Count;

            // calculating the worst fitness value in the current generation (higher fitness value is worse)
            currentGeneration_WorstFitness = int.MinValue;
            for (int i = 0; i < generationPopulation; i++)
            {
                if (currentGeneration_WorstFitness < currentGeneration.chromosomes[i].fitness)
                    currentGeneration_WorstFitness = currentGeneration.chromosomes[i].fitness;
            }

            alpha = alphaValue;
            beta = betaValue;
            double sumOfFitnessValues = 0.0;
            for (int i = 0; i < generationPopulation; i++)
            {
                sumOfFitnessValues += 1.0 / currentGeneration.chromosomes[i].fitness;
                if (currentGeneration.chromosomes[i].fitness <= currentGeneration_WorstFitness / alpha)
                    sumOfFitnessValues += beta * (1.0 / currentGeneration.chromosomes[i].fitness);
            }

            portion = 0;
            portionPercentage = new double[generationPopulation];
            for (int i = 0; i < generationPopulation; i++)
            {
                portion += (1.0 / currentGeneration.chromosomes[i].fitness) / sumOfFitnessValues; // chance of each chromosome 
                // if a chromosome has a better fitness than the worst one "by equal or more than beta times"
                if (currentGeneration.chromosomes[i].fitness <= currentGeneration_WorstFitness / alpha)
                    // give it a higher chance of selection "by alpha times"
                    portion += beta * ((1.0 / currentGeneration.chromosomes[i].fitness) / sumOfFitnessValues);

                portionPercentage[i] = portion;
            }

            winnerChromosomes = new List<Chromosome>();
            for (int i = 0; i < generationPopulation; i++) // selecting winner chromosomes based on their fitness value, randomly 
            {
                double probability = rnd.NextDouble();

                for (int j = 0; j < generationPopulation; j++)
                {
                    if (probability <= portionPercentage[j])
                    {
                        winnerChromosomes.Add(currentGeneration.chromosomes[j]);
                        break;
                    }
                }
            }

            return new Generation(winnerChromosomes);

        }
    }

    class TournamentSelection
    {
        public static int generationPopulation;
        public static List<Chromosome> Tournament;
        public static List<Chromosome> winnerChromosomes;


        public static int randomNumber;
        public static Random rnd = new Random();

        public static Generation ApplyTournamentSelection(Generation currentGeneration, int k)
        {
            generationPopulation = currentGeneration.chromosomes.Count;
            winnerChromosomes = new List<Chromosome>();
            for (int i = 0; i < generationPopulation; i++) // selecting single new chromosomes (by winning the tournament), as many as number of the population
            {
                Tournament = new List<Chromosome>();
                for (int j = 0; j < k; j++) // selecting k random chromosomes
                {
                    randomNumber = rnd.Next(generationPopulation); // choosing the index of chromosome randomly
                    Tournament.Add(currentGeneration.chromosomes[randomNumber]); // adding a random chromosome from last generation to tournament
                }
                Tournament = Tournament.OrderBy(chromosome => chromosome.fitness).ToList(); // ordering the chromosomes in the tournament based on their fitness value

                winnerChromosomes.Add(Tournament[0]); // creating the list of winner chromosomes 
            }

            return new Generation(winnerChromosomes);
        }
    }

    class Crossover
    {
        public static Chromosome[] ApplyCrossover(Chromosome parentChromosome1, Chromosome parentChromosome2, string crossoverMethod)
        {
            if (crossoverMethod == "singlePoint")
            {
                if(parentChromosome1.creationProcess == "standard")
                    return SinglePointCrossover.ApplySinglePointCrossover(parentChromosome1, parentChromosome2);
                else if(parentChromosome1.creationProcess == "permutation")
                    return Repair.ApplyRepair(
                        SinglePointCrossover.ApplySinglePointCrossover(parentChromosome1, parentChromosome2));
            }
            else if (crossoverMethod == "uniform")
            {
                if (parentChromosome1.creationProcess == "standard")
                    return UniformCrossover.ApplyUniformCrossover(parentChromosome1, parentChromosome2);
                else if (parentChromosome1.creationProcess == "permutation")
                    return Repair.ApplyRepair(
                        UniformCrossover.ApplyUniformCrossover(parentChromosome1, parentChromosome2));
            }
            else if (crossoverMethod == "order1")
            {
                if (parentChromosome1.creationProcess == "standard")
                    return Order1Crossover.ApplyOrder1Crossover(parentChromosome1, parentChromosome2); // it's not neccessary to have order1 on standard, but it's there btw
                else if (parentChromosome1.creationProcess == "permutation")
                    return Order1Crossover.ApplyOrder1Crossover(parentChromosome1, parentChromosome2); // order1 on permutation doesn't need repair mechanism
            }
            
            return null; 
        }
    }

    class SinglePointCrossover 
    {
        // input chromosomes
        public static int[] parent1;
        public static int[] parent2;

        // output chromosomes with crossover applied
        public static int[] child1;
        public static int[] child2;

        public static int crossoverPoint;

        public static Chromosome[] newChromosomes;

        public static Random rnd = new Random();

        public static Chromosome[] ApplySinglePointCrossover(Chromosome inputChromosome1, Chromosome inputChromosome2)
        {
            parent1 = (int[])inputChromosome1.chromosome.Clone();
            parent2 = (int[])inputChromosome2.chromosome.Clone();

            int length = parent1.Length;

            // initialize
            child1 = new int[length];
            child2 = new int[length];

            crossoverPoint = rnd.Next(length); // make a crossover point

            for (int i = 0; i < length; ++i) // making two childs 
            {
                // left side crossover
                if (i < crossoverPoint) 
                {
                    child1[i] = parent1[i];
                    child2[i] = parent2[i];
                }
                // right side crossover
                else
                {
                    child1[i] = parent2[i];
                    child2[i] = parent1[i];
                }
            }

            // creating 2 new child chromosomes and sending them back to genetic algorithm
            return newChromosomes = new Chromosome[] {
                new Chromosome(child1, inputChromosome1.creationProcess),
                new Chromosome(child2, inputChromosome1.creationProcess)};
        }
    }

    class UniformCrossover
    {
        // input chromosomes
        public static int[] parent1;
        public static int[] parent2;

        // output chromosomes with crossover applied
        public static int[] child1;
        public static int[] child2;

        public static double crossoverProbability = 0.5; // default: 0.5

        public static Chromosome[] newChromosomes;

        public static Random rnd = new Random();

        public static Chromosome[] ApplyUniformCrossover(Chromosome inputChromosome1, Chromosome inputChromosome2)
        {
            parent1 = (int[])inputChromosome1.chromosome.Clone();
            parent2 = (int[])inputChromosome2.chromosome.Clone();

            int length = parent1.Length;

            // initialize
            child1 = new int[length];
            child2 = new int[length];

            for (int i = 0; i < length; i++) // making 2 childs
            {
                // applying crossover in child chromosomes
                if (rnd.NextDouble() < crossoverProbability) 
                {
                    child1[i] = parent2[i];
                    child2[i] = parent1[i];
                }
                // transferring genes from parent to child without modification 
                else
                {
                    child1[i] = parent1[i];
                    child2[i] = parent2[i];
                }
            }

            // creating 2 new child chromosomes and sending them back to genetic algorithm
            return newChromosomes = new Chromosome[] {
                new Chromosome(child1, inputChromosome1.creationProcess),
                new Chromosome(child2, inputChromosome1.creationProcess)};
        }
    }

    class Order1Crossover
    {
        // input chromosomes
        public static int[] parent1;
        public static int[] parent2;

        // output chromosomes with crossover applied
        public static int[] child;

        // output chromosomes within a single function 
        public static Chromosome child1;
        public static Chromosome child2;

        public static Chromosome[] newChromosomes;

        public static Random rnd = new Random();

        public static Chromosome[] ApplyOrder1Crossover(Chromosome inputChromosome1, Chromosome inputChromosome2)
        {
            child1 = Order1Crossover_SingleOutput(inputChromosome1, inputChromosome2);
            child2 = Order1Crossover_SingleOutput(inputChromosome2, inputChromosome1);

            // creating 2 new child chromosomes and sending them back to genetic algorithm
            return newChromosomes = new Chromosome[] {child1, child2};
        }

        public static Chromosome Order1Crossover_SingleOutput(Chromosome inputChromosome1, Chromosome inputChromosome2)
        {
            parent1 = (int[])inputChromosome1.chromosome.Clone(); 
            parent2 = (int[])inputChromosome2.chromosome.Clone();

            int length = parent1.Length;

            // get 2 random ints between 0 and size of array
            int randomNum1 = rnd.Next(length);
            int randomNum2 = rnd.Next(length);

            while (randomNum1 >= randomNum2) // to make sure that randomNum1 < randomNum2
            {
                randomNum1 = rnd.Next(length);
                randomNum2 = rnd.Next(length);
            }

            // child chromosome creation with parent chromosomes length
            child = new int[length];
            for (int i = 0; i < child.Length; i++) // initial elements with -1 values
            {
                child[i] = -1;
            }

            // copy elements between randomNum1, randomNum2 from parent1 into child
            for (int i = randomNum1; i <= randomNum2; i++)
            {
                child[i] = parent1[i];
            }

            // array to hold elements of parent1 which are not in child yet
            int[] y = new int[length - (randomNum2 - randomNum1) - 1];
            int j = 0;
            for (int i = 0; i < length; i++)
            {
                if (!arrayContains(child, parent1[i]))
                {
                    y[j] = parent1[i];
                    j++;
                }
            }

            // rotate parent2
            // number of places is the same as the number of elements after randomNum2
            int[] copy = (int[])parent2.Clone();
            RotateToRight(copy, length - randomNum2 - 1);

            // now order the elements in y according to their order in parent2 
            int[] y1 = new int[length - (randomNum2 - randomNum1) - 1];
            j = 0;
            for (int i = 0; i < length; i++)
            {
                if (arrayContains(y, copy[i]))
                {
                    y1[j] = copy[i];
                    j++;
                }
            }

            // now copy the remaining elements (i.e. remaining in parent1) into child 
            // according to their order in parent2 .. starting after randomNum2
            j = 0;
            for (int i = 0; i < y1.Length; i++)
            {
                int ci = (randomNum2 + i + 1) % length; // current index
                child[ci] = y1[i];
            }

            return new Chromosome(child, inputChromosome1.creationProcess);
        }

        // checks if an int[] array contains a specific int element
        public static bool arrayContains(int[] arr, int e)
        {
            bool result = false;

            for (int i = 0; i < arr.Length; i++)
            {
                if (arr[i] == e)
                    result = true;
            }
            return result;
        }

        // rotating the elements after randomNum2 
        public static void RotateToRight(int[] input, int numberOfShifts) // function to shift an array to right, amount is specified with numberOfShifts
        {
            for (var i = 0; i < numberOfShifts; i += 1)
            {
                RotateRightOne(input);
            }
        }

        public static void RotateRightOne(int[] input) // function to shift an array to right, only one time
        {
            var last = input.Length - 1;
            for (var i = 0; i < last; i += 1)
            {
                input[i] ^= input[last];
                input[last] ^= input[i];
                input[i] ^= input[last];
            }
        }

    }

    class Repair
    {
        public static int[] input; // input chromosome
        public static Chromosome outputChromosome; // output chromosome with repair modification applied

        // a template chromosome with all possible gene values and no duplication || i.e. for N:3 => {0, 1, 2}
        public static int[] template;

        // stores how many times a specefic allele is repeated in a chromosome {Key: allele, Value: number of repetitions}
        public static Dictionary<int, int> duplicationCounter = new Dictionary<int, int>();
        // stores a list of duplicated values, each repeated value only exists once in this list
        public static List<int> duplicates = new List<int>();
        // for using in input chromosome reconstruction, with this, the first occurence of a repeated value won't change to -1
        public static Dictionary<int, int> duplicateCountedValue = new Dictionary<int, int>();
        // a list of values that are not present in the input chromosme 
        public static List<int> remainingValues = new List<int>();
        // for using in input chromosome reconstruction, with this, each value from remainingValues, will be used only once
        public static List<int> usedRemainingValues = new List<int>();

        public static Chromosome child1;
        public static Chromosome child2;
        public static Chromosome[] newChromosomes;

        public static Random rnd = new Random();

        public static Chromosome[] ApplyRepair(Chromosome[] inputChromosome)
        {
            child1 = Repair_SingleOutput(inputChromosome[0]);
            child2 = Repair_SingleOutput(inputChromosome[1]);

            // creating 2 new child chromosomes and sending them back to genetic algorithm
            return newChromosomes = new Chromosome[] { child1, child2 };
        }

        public static Chromosome Repair_SingleOutput(Chromosome inputChromosome)
        {
            input = (int[])inputChromosome.chromosome.Clone();

            int length = input.Length;

            template = TemplateCreation(length);

            usedRemainingValues = new List<int>(); // creating a new usedRemainingValues list, each time we use the function

            foreach (int allele in input) // finding duplicates in the input chromosome
            {
                if (!duplicationCounter.ContainsKey(allele))
                    duplicationCounter.Add(allele, 0);
                duplicationCounter[allele]++;
            }
            // add each value to the list of duplicates if it appears 2 times or more in the chromosome
            duplicates = duplicationCounter.Where(c => c.Value >= 2).Select(c => c.Key).ToList();

            // traversing the input chromosome and changing the index value of duplicates to -1 
            for (int i = 0; i < input.Length; i++)
            {
                if (duplicates.Contains(input[i]))
                {
                    // with this, the first occurence of a repeated value won't change to -1, only the second, third.... occurrences will change 
                    // this is to preserve the order of appearance of duplicate values 
                    if (!duplicateCountedValue.ContainsKey(input[i]))
                        duplicateCountedValue.Add(input[i], 0);
                    duplicateCountedValue[input[i]]++;

                    if (duplicateCountedValue[input[i]] >= 2)
                        input[i] = -1;
                }
            }

            // making a remaining values list which are not present in the input chromosme, but are in the template chromosme
            remainingValues = template.Except(input).ToList();

            // adding the remaining values to the input chromosome
            for (int i = 0; i < input.Length; i++)
            {
                if (input[i] == -1) // if a duplicate index found 
                {
                    for (int j = 0; j < remainingValues.Count; j++) // look for all the remaining values
                    {
                        if (!usedRemainingValues.Contains(remainingValues[j])) // if a remaining value never used, use it
                        {
                            input[i] = remainingValues[j];
                            usedRemainingValues.Add(remainingValues[j]); // add the remaining value to the used ones 
                            break;
                        }
                    }
                }
            }

            return outputChromosome = new Chromosome(input, inputChromosome.creationProcess);
        }

        // a function to create a template chromosome with all possible gene values and no duplication || i.e. for N:3 => {0, 1, 2}
        public static int[] TemplateCreation(int length)
        {
            int[] temp = new int[length];
            for (int i = 0; i < length; i++)
                temp[i] = i;

            return temp;
        }
    }

    class Mutation
    {
        public static int[] outputChromosome;

        public static Random rnd = new Random();

        public static Chromosome ApplyMutation(Chromosome inputChromosome)
        {
            if (inputChromosome.creationProcess == "standard")
                return StandardMutation(inputChromosome);
            else 
                return SwapMutation(inputChromosome);
        }

        public static Chromosome StandardMutation(Chromosome inputChromosome)
        {
            // Gives a random number between 0,...,N-1 (we use this number to determine which queen's position will be changed)
            int randomQueen = rnd.Next(inputChromosome.chromosome.Length);
            // Gives a random number between 0,...,N-1 (we use this number to change queen's position horizontally)
            int randomHorizontalPosition = rnd.Next(inputChromosome.chromosome.Length);

            outputChromosome = (int[])inputChromosome.chromosome.Clone();
            // Only one queen's position is changed horizontally (repetetive satates may occure)
            outputChromosome[randomQueen] = randomHorizontalPosition;

            return new Chromosome(outputChromosome, inputChromosome.creationProcess);
        }

        public static Chromosome SwapMutation(Chromosome inputChromosome)
        {
            int randomRow1 = rnd.Next(inputChromosome.chromosome.Length); // Gives a random number between 0,...,N-1 
            int randomRow2 = rnd.Next(inputChromosome.chromosome.Length); // Gives a random number between 0,...,N-1 

            outputChromosome = (int[])inputChromosome.chromosome.Clone();
            // Only one queen's position is changed horizontally (repetetive satates will not occure)
            outputChromosome[randomRow1] = inputChromosome.chromosome[randomRow2];
            outputChromosome[randomRow2] = inputChromosome.chromosome[randomRow1];

            return new Chromosome(outputChromosome, inputChromosome.creationProcess);
        }
    }

    class Fitness
    {
        public static int numberOfAttacks; // number of attacks, which can be in the form of vertical and diagonal attacks between queens (Goal state has 0 numberOfAttacks)
        public static int FitnessFunction(int[] chromosome) // calculates the number of attacks and returns the number (the less, the better)
        {
            numberOfAttacks = 0; // resetting the value 

            for (int i = 0; i < chromosome.Length; i++)
            {
                for (int j = i + 1; j < chromosome.Length; j++)
                {
                    if (chromosome[i] == chromosome[j]) // calculating the number of vertical attacks
                        numberOfAttacks++;
                    if (Math.Abs(i - j) == Math.Abs(chromosome[i] - chromosome[j])) // calculating the number of diagonal attacks
                        numberOfAttacks++;
                }
            }
            return numberOfAttacks;
        }
    }
}
