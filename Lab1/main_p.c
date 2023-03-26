#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

#define DT 0.05

typedef struct
{
    double x, y;
} vector;

typedef struct
{
    vector first;
    vector second;
} VectorPair;

typedef struct Task {
    int i;
} Task;

typedef struct PairTask {
    int i, j;
} PairTask;


int bodies, timeSteps, pair_num, threads_num = 4;
double *masses, GravConstant;
double *grav_mass_constants;
vector *positions, *velocities, *accelerations;

Task* singleTaskQueue;
int singleTaskCount = 0, singleLen = 0;
PairTask* pairTaskQueue;
int pairTaskCount = 0, pairLen = 0;

int cycle_computed = 0, cycle_computed_single = 0;

// thread pool
pthread_mutex_t mutexQueue, mutexSingleQueue;
pthread_cond_t condQueue, condSingleQueue;

// mutex per body
pthread_mutex_t* mutexes;

pthread_t* threads_pool;
pthread_t* threads_pool_2;

void submitTask(Task task) {
    // pthread_mutex_lock(&mutexQueue);
    singleTaskQueue[singleLen] = task;
    singleLen++;
    // pthread_mutex_unlock(&mutexQueue);
    // pthread_cond_signal(&condQueue);
}

void resetTaskPool() {
    pthread_mutex_lock(&mutexSingleQueue);
    singleTaskCount = bodies;
    cycle_computed_single = 0;
    pthread_mutex_unlock(&mutexSingleQueue);
    pthread_cond_signal(&condSingleQueue);
}

void submitPairTask(PairTask task) {
    // pthread_mutex_lock(&mutexQueue);
    pairTaskQueue[pairLen] = task;
    pairLen++;
    // pthread_mutex_unlock(&mutexQueue);
    // pthread_cond_signal(&condQueue);
}

void resetPairTaskPool() {
    pthread_mutex_lock(&mutexQueue);
    pairTaskCount = pair_num;
    cycle_computed = 0;
    pthread_mutex_unlock(&mutexQueue);
    pthread_cond_signal(&condQueue);
}

void* startThread(void* func) {
    // printf("SINGLE TASK COUNT %d\n", singleTaskCount);
    while (1) {
        Task task;

        pthread_mutex_lock(&mutexSingleQueue);
        while (singleTaskCount == 0) {
            cycle_computed_single = 1;
            pthread_cond_wait(&condSingleQueue, &mutexSingleQueue);
        }

        task = singleTaskQueue[singleTaskCount - 1];
        // ("SINGLE TASK BEFORE MINUS %d\n", singleTaskCount);
        singleTaskCount--;
        // printf("SINGLE COUNT AFTER MINUS %d", singleTaskCount);
        pthread_mutex_unlock(&mutexSingleQueue);
        void (*fun_ptr)(int) = func;
        (*fun_ptr)(task.i);
    }
    return NULL;
}

void* startPairThread(void* func) {
    while (1) {
        PairTask task;

        pthread_mutex_lock(&mutexQueue);
        while (pairTaskCount == 0) {
            cycle_computed = 1;
            pthread_cond_wait(&condQueue, &mutexQueue);
        }

        task = pairTaskQueue[pairTaskCount - 1];
        //printf("pairTaskCount=%d\n", pairTaskCount);
        //printf("taking Task[%d]{%d, %d}\n", pairTaskCount - 1, task.i, task.j);
        pairTaskCount--;
        //printf("pairTaskCount=%d\n", pairTaskCount);
        pthread_mutex_unlock(&mutexQueue);
        void (*fun_ptr)(int, int) = func;
        (*fun_ptr)(task.i, task.j);
    }
    return NULL;
}

vector addVectors(vector a, vector b)
{
    vector c = {a.x + b.x, a.y + b.y};
    return c;
}

vector scaleVector(double b, vector a)
{
    vector c = {b * a.x, b * a.y};
    return c;
}

vector subtractVectors(vector a, vector b)
{
    vector c = {a.x - b.x, a.y - b.y};
    return c;
}

double mod(vector a)
{
    return sqrt(a.x * a.x + a.y * a.y);
}

int computePairsNum()
{
    if (bodies < 2) return 0;
    int numerator = 1;
    int i;
    for (i = bodies - 1; i <= bodies; ++i)
        numerator *= i;
    return numerator / 2;
}

void computeGravityMassConstants()
{
    for (int i = 0; i < bodies; ++i) {
        grav_mass_constants[i] = masses[i] * GravConstant;
    }
}

void UpdatePair(int i, int j);
void UpdateSingle(int i);
void computeVelocity(int i);
void computePosition(int i);



void initiateSystem(char *fileName)
{
    int i, j;
    FILE *fp = fopen(fileName, "r");

    fscanf(fp, "%lf%d%d", &GravConstant, &bodies, &timeSteps);

    masses = (double *)malloc(bodies * sizeof(double));
    grav_mass_constants = (double *)malloc(bodies * sizeof(double));
    positions = (vector *)malloc(bodies * sizeof(vector));
    velocities = (vector *)malloc(bodies * sizeof(vector));
    accelerations = (vector *)malloc(bodies * sizeof(vector));

    mutexes = (pthread_mutex_t *)malloc(bodies * sizeof(pthread_mutex_t));
    for (i = 0; i < bodies; ++i)
        pthread_mutex_init(&(mutexes[i]), NULL);
    pthread_mutex_init(&mutexQueue, NULL);
    pthread_cond_init(&condQueue, NULL);
    pthread_mutex_init(&mutexSingleQueue, NULL);
    pthread_cond_init(&condSingleQueue, NULL);

    // pools of tasks
    pair_num = computePairsNum();
    singleTaskQueue = (Task *)malloc(bodies * sizeof(Task));
    pairTaskQueue = (PairTask *)malloc(pair_num * sizeof(PairTask));
    // fill pool of tasks
    pairTaskCount = 0;
    //printf("create task pool\n");

    for (i = 0; i < bodies; i++) {
        submitTask((Task) {i});
    }

    for (i = 0; i < bodies - 1; i++)
        for (j = i + 1; j < bodies; j++) {
            submitPairTask((PairTask){i , j});
        }
    // pool of threads
    threads_pool = (pthread_t *)malloc(threads_num * sizeof(pthread_t));
    threads_pool_2 = (pthread_t *)malloc(threads_num * sizeof(pthread_t));
    //printf("create threads\n");
    for (i = 0; i < threads_num; i++) {
        void (*fun_ptr)(int, int) = &UpdatePair;
        if (pthread_create(&threads_pool[i], NULL, &startPairThread, fun_ptr) != 0) {
            perror("Failed to create the thread");
        }
    }

    for (i = 0; i < threads_num; i++) {
        void (*fun_ptr)(int) = &UpdateSingle;
        if (pthread_create(&threads_pool_2[i], NULL, &startThread, fun_ptr) != 0) {
            perror("Failed to create the thread");
        }
    }

    // Take params from file
    for (i = 0; i < bodies; i++)
    {
        fscanf(fp, "%lf", &masses[i]);
        fscanf(fp, "%lf%lf", &positions[i].x, &positions[i].y);
        fscanf(fp, "%lf%lf", &velocities[i].x, &velocities[i].y);
    }

    computeGravityMassConstants();

    fclose(fp);
}

void destructSystem()
{
    int i;
    for (i = 0; i < bodies; ++i)
        pthread_mutex_destroy(&(mutexes[i]));

    pthread_mutex_destroy(&mutexQueue);
    pthread_cond_destroy(&condQueue);
    pthread_mutex_destroy(&mutexSingleQueue);
    pthread_cond_destroy(&condSingleQueue);

    for (i = 0; i < threads_num; ++i) {
        pthread_cancel(threads_pool[i]);
        pthread_cancel(threads_pool_2[i]);
    }
}

void resolveCollisions()
{
    int i, j;

    for (i = 0; i < bodies - 1; i++)
    {
        for (j = i + 1; j < bodies; j++)
        {
            if (positions[i].x == positions[j].x && positions[i].y == positions[j].y)
            {
                vector temp = velocities[i];
                velocities[i] = velocities[j];
                velocities[j] = temp;
            }
        }
    }
}

// ReadOnly
vector computeCommonPart(int i, int j) {
    vector sub = subtractVectors(positions[j], positions[i]);
    double len = mod(sub);
    double k = 1 / pow(len, 3);
    return scaleVector(k, sub);
}

// ReadOnly
VectorPair ComputeIndependentPair(int i, int j)
{
    vector common_part = computeCommonPart(i, j);
    return (VectorPair) {
            scaleVector( 1 * grav_mass_constants[j], common_part),
            scaleVector(-1 * grav_mass_constants[i], common_part)
    };
}


void UpdatePair(int i, int j)
{
    //printf("running Task{%d, %d}\n", i, j);
    VectorPair independent_pair = ComputeIndependentPair(i, j);
    // vector add_part_ij = independent_pair.first;
    // vector add_part_ji = independent_pair.second;
    // lock for i'th body

    pthread_mutex_lock(&(mutexes[i]));

    //printf("accelerations[%d] = {%lf, %lf}\n", i, accelerations[i].x, accelerations[i].y);
    //printf("independent_pair.i[%d] = {%lf, %lf}\n", i, add_part_ij.x, add_part_ij.y);
    accelerations[i] = addVectors(accelerations[i], independent_pair.first);

    pthread_mutex_unlock(&(mutexes[i]));

    // lock for j'th body

    pthread_mutex_lock(&(mutexes[j]));

    //("accelerations[%d] = {%lf, %lf}\n", j, accelerations[j].x, accelerations[j].y);
    //printf("independent_pair.j[%d] = {%lf, %lf}\n", j, add_part_ji.x, add_part_ji.y);
    accelerations[j] = addVectors(accelerations[j], independent_pair.second);

    pthread_mutex_unlock(&(mutexes[j]));
}


void UpdateSingle(int i)
{
    computePosition(i);
    computeVelocity(i);
}

void resetAccelerations()
{
    for (int i = 0; i < bodies; i++)
    {
        accelerations[i].x = 0.f;
        accelerations[i].y = 0.f;
    }
}

void computeAccelerationsInParallel()
{
    // create task pool
    // int i, j;
    resetPairTaskPool();
    // join threads
    //printf("join threads\n");
    while(!cycle_computed) {}
}

void computeAccelerations()
{
    int i, j;

    for (i = 0; i < bodies; i++)
    {
        accelerations[i].x = 0;
        accelerations[i].y = 0;
        for (j = 0; j < bodies; j++)
        {
            if (i != j)
            {
                vector add_part =
                        scaleVector(
                                GravConstant * masses[j] /
                                pow(mod(subtractVectors(positions[i], positions[j])), 3),
                                subtractVectors(positions[j], positions[i])
                        );
                accelerations[i] = addVectors(
                        accelerations[i],
                        add_part
                );
            }
        }
    }
}

void computeVelocity(int i)
{
    // printf("Compute velocity\n");
    velocities[i] = addVectors(velocities[i], scaleVector(DT, accelerations[i]));
}

void computePosition(int i)
{
    // printf("Compute position\n");
    positions[i] = addVectors(positions[i], scaleVector(DT, velocities[i]));
}

void computeVelocities()
{
    for(int i = 0; i < bodies; ++i)
        computeVelocity(i);
}

void computePositions()
{
    for(int i = 0; i < bodies; ++i)
        computePosition(i);
}


void computeVelocitiesAndPositionsInParallel()
{
    resetTaskPool();
    while(!cycle_computed_single) {}
}


void simulate()
{
    // computeAccelerations();
    resetAccelerations();
    computeAccelerationsInParallel();
    computeVelocitiesAndPositionsInParallel();
//    computePositions();  // для каждого тела
//    computeVelocities(); // для каждого тела
    resolveCollisions(); // не параллельно
}

int main(int argC, char *argV[])
{
    int i, j;
    if (argC != 2)
    {
        printf("Usage : %s <file name containing system configuration data> <threads_num>", argV[0]);
        return 0;
    }

    initiateSystem(argV[1]);

    printf("Body   :     x              y           vx              vy   ");
    double time_spent = 0.0;
    clock_t begin = clock();
    for (i = 0; i < timeSteps; i++)
    {
        printf("\nCycle %d\n", i + 1);
        //double time_spent_one_it = 0.0;
        //clock_t begin_one_it = clock();
        simulate();
        //clock_t end_one_it = clock();
        //time_spent_one_it += (double)(end_one_it - begin_one_it) / CLOCKS_PER_SEC;
        for (j = 0; j < bodies; j++)
            printf("Body %d : %lf\t%lf\t%lf\t%lf\n", j + 1, positions[j].x, positions[j].y, velocities[j].x, velocities[j].y);
        //printf("\nThe elapsed time is %.10lf seconds\n", time_spent_one_it);
    }
    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The elapsed time is %f seconds", time_spent);
    destructSystem();

    return 0;
}