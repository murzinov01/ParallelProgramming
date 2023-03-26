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


int bodies, timeSteps, pair_num, threads_num = 8;
double *masses, GravConstant;
double *grav_mass_constants;
vector *positions, *velocities, *accelerations;

Task* singleTaskQueue;
int singleTaskCount = 0;
PairTask* pairTaskQueue;
int pairTaskCount = 0;

// thread pool
pthread_mutex_t mutexPairQueue, mutexSingleQueue;
pthread_cond_t condPairQueue, condSingleQueue;

pthread_mutex_t mutexCycle;
pthread_cond_t condCycle;

// mutex per body
pthread_mutex_t* mutexes;

pthread_t* threads_pair_pool;
pthread_t* threads_single_pool;

void resetTaskPool() {
    pthread_mutex_lock(&mutexSingleQueue);
    singleTaskCount = bodies;
    pthread_mutex_unlock(&mutexSingleQueue);
    pthread_cond_signal(&condSingleQueue);
}

void resetPairTaskPool() {
    pthread_mutex_lock(&mutexPairQueue);
    pairTaskCount = pair_num;
    pthread_mutex_unlock(&mutexPairQueue);
    pthread_cond_signal(&condPairQueue);
}

void* startThread(void* func) {
    int ran = 0;
    while (1) {
        Task task;

        pthread_mutex_lock(&mutexSingleQueue);
        while (singleTaskCount == 0) {
            if (ran) {
                pthread_cond_signal(&condCycle);
            }
            pthread_cond_wait(&condSingleQueue, &mutexSingleQueue);
        }

        task = singleTaskQueue[singleTaskCount - 1];
        singleTaskCount--;
        pthread_mutex_unlock(&mutexSingleQueue);
        void (*fun_ptr)(int) = func;
        (*fun_ptr)(task.i);
        ran = 1;
    }
    return NULL;
}

void* startPairThread(void* func) {
    int ran = 0;
    while (1) {
        PairTask task;

        pthread_mutex_lock(&mutexPairQueue);
        while (pairTaskCount == 0) {
            if (ran) {
                pthread_cond_signal(&condCycle);
            }
            pthread_cond_wait(&condPairQueue, &mutexPairQueue);
        }

        task = pairTaskQueue[pairTaskCount - 1];
        pairTaskCount--;
        pthread_mutex_unlock(&mutexPairQueue);
        void (*fun_ptr)(int, int) = func;
        (*fun_ptr)(task.i, task.j);
        ran = 1;
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

void FillTasks() {
    int i, j;
    for (i = 0; i < bodies; i++)
        singleTaskQueue[i] = (Task) {i};

    size_t pair_len = 0;
    for (i = 0; i < bodies - 1; i++)
        for (j = i + 1; j < bodies; j++)
            pairTaskQueue[pair_len++] = (PairTask){i , j};
}

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
    pthread_mutex_init(&mutexPairQueue, NULL);
    pthread_cond_init(&condPairQueue, NULL);
    pthread_mutex_init(&mutexSingleQueue, NULL);
    pthread_cond_init(&condSingleQueue, NULL);

    // pools of tasks
    pair_num = computePairsNum();
    singleTaskQueue = (Task *)malloc(bodies * sizeof(Task));
    pairTaskQueue = (PairTask *)malloc(pair_num * sizeof(PairTask));
    // fill pool of tasks
    FillTasks();

    // pool of threads
    threads_pair_pool = (pthread_t *)malloc(threads_num * sizeof(pthread_t));
    threads_single_pool = (pthread_t *)malloc(threads_num * sizeof(pthread_t));

    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS,NULL);
    for (i = 0; i < threads_num; i++) {
        void (*fun_ptr)(int, int) = &UpdatePair;
        if (pthread_create(&threads_pair_pool[i], NULL, &startPairThread, fun_ptr) != 0) {
            perror("Failed to create the thread");
        }
    }

    for (i = 0; i < threads_num; i++) {
        void (*fun_ptr)(int) = &UpdateSingle;
        if (pthread_create(&threads_single_pool[i], NULL, &startThread, fun_ptr) != 0) {
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

    for (i = 0; i < threads_num; ++i) {
        pthread_cancel(threads_pair_pool[i]);
        pthread_cancel(threads_single_pool[i]);
    }

    pthread_mutex_destroy(&mutexPairQueue);
    pthread_cond_destroy(&condPairQueue);
    pthread_mutex_destroy(&mutexSingleQueue);
    pthread_cond_destroy(&condSingleQueue);
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
    VectorPair independent_pair = ComputeIndependentPair(i, j);
    // lock for i'th body
    pthread_mutex_lock(&(mutexes[i]));
    accelerations[i] = addVectors(accelerations[i], independent_pair.first);
    pthread_mutex_unlock(&(mutexes[i]));

    // lock for j'th body
    pthread_mutex_lock(&(mutexes[j]));
    accelerations[j] = addVectors(accelerations[j], independent_pair.second);
    pthread_mutex_unlock(&(mutexes[j]));
}

void computeAccelerationsInParallel()
{
    pthread_mutex_lock(&mutexCycle);
    resetPairTaskPool();
    pthread_cond_wait(&condCycle, &mutexCycle);
    pthread_mutex_unlock(&mutexCycle);
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

void resetAccelerations()
{
    for (int i = 0; i < bodies; i++)
    {
        accelerations[i].x = 0.f;
        accelerations[i].y = 0.f;
    }
}

void computeVelocity(int i)
{
    velocities[i] = addVectors(velocities[i], scaleVector(DT, accelerations[i]));
}

void computePosition(int i)
{
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

void UpdateSingle(int i)
{
    computePosition(i);
    computeVelocity(i);
}

void computeVelocitiesAndPositionsInParallel()
{
    pthread_mutex_lock(&mutexCycle);
    resetTaskPool();
    pthread_cond_wait(&condCycle, &mutexCycle);
    pthread_mutex_unlock(&mutexCycle);
}

void simulateInParallel()
{
    resetAccelerations();
    computeAccelerationsInParallel();
    computeVelocitiesAndPositionsInParallel();
    resolveCollisions(); // не параллельно
}

void simulate()
{
    computeAccelerations();
    computePositions();
    computeVelocities();
    resolveCollisions();
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

        simulate();
        // simulateInParallel();

        for (j = 0; j < bodies; j++)
            printf("Body %d : %lf\t%lf\t%lf\t%lf\n", j + 1, positions[j].x, positions[j].y, velocities[j].x, velocities[j].y);
    }
    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The elapsed time is %f seconds", time_spent);

    destructSystem();

    return 0;
}