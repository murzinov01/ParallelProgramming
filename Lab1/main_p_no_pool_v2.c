#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>

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


int bodies, timeSteps, pair_num, threads_num = 24;
double *masses, GravConstant;
double *grav_mass_constants;
vector *positions, *velocities, *accelerations;

Task* singleTaskQueue;
int singleTaskCount = 0, singleLen = 0;
PairTask* pairTaskQueue;
int pairTaskCount = 0, pairLen = 0;

int** tasksForThread;
size_t* tasksForThreadLens;

// mutex per body
pthread_mutex_t* mutexes;
pthread_t* threads_pool;

pthread_mutex_t end_mutex;

void submitTask(Task task) {
    singleTaskQueue[singleTaskCount] = task;
    singleTaskCount++;
}

void submitPairTask(PairTask task) {
    pairTaskQueue[pairLen] = task;
    pairLen++;
}

void UpdatePair(int, int);

void* runPairThread(void* thread_id) {
    int thread_id_ = *((int *) thread_id);
    for (int i = 0; i < tasksForThreadLens[thread_id_]; ++i) {
        int task_id = tasksForThread[thread_id_][i];
        PairTask task = pairTaskQueue[task_id];
        UpdatePair(task.i, task.j);
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

void initiateSystem(char *fileName)
{
    int i, j;
    FILE *fp = fopen(fileName, "r");

    fscanf(fp, "%lf%d%d", &GravConstant, &bodies, &timeSteps);

    pair_num = computePairsNum();

    masses = (double *)malloc(bodies * sizeof(double));
    grav_mass_constants = (double *)malloc(bodies * sizeof(double));
    positions = (vector *)malloc(bodies * sizeof(vector));
    velocities = (vector *)malloc(bodies * sizeof(vector));
    accelerations = (vector *)malloc(bodies * sizeof(vector));

    singleTaskQueue = (Task *)malloc(bodies * sizeof(Task));
    pairTaskQueue = (PairTask *)malloc(pair_num * sizeof(PairTask));

    pthread_mutex_init(&end_mutex, NULL);
    mutexes = (pthread_mutex_t *)malloc(bodies * sizeof(pthread_mutex_t));
    for (i = 0; i < bodies; ++i)
        pthread_mutex_init(&(mutexes[i]), NULL);

    tasksForThreadLens = (size_t *)malloc(threads_num * sizeof(size_t));
    tasksForThread = (int **)malloc(threads_num * sizeof(int*));
    int cur_thread = 0;
    for (i = 0; i < threads_num; ++i)
        tasksForThreadLens[i] = 0;
    for (i = 0; i < pair_num; ++i)
        ++tasksForThreadLens[cur_thread++ % threads_num];
    for (i = 0; i < threads_num; ++i)
        tasksForThread[i] = (int *)malloc(tasksForThreadLens[i] * sizeof(int));
    int cur_task = 0;
    for (i = 0; i < threads_num; ++i) {
        for (j = 0; j < tasksForThreadLens[i]; ++j) {
            tasksForThread[i][j] = cur_task++;
        }
    }

    // fill pool of tasks
    pairTaskCount = 0;
    printf("create task pool\n");
    for (i = 0; i < bodies - 1; i++)
        for (j = i + 1; j < bodies; j++) {
            submitPairTask((PairTask){i , j});
        }
    // pool of threads
    threads_pool = (pthread_t *)malloc(threads_num * sizeof(pthread_t));
    printf("create threads\n");

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

    for (i = 0; i < threads_num; ++i)
        pthread_cancel(threads_pool[i]);

    pthread_mutex_destroy(&end_mutex);
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
    // printf("running Task{%d, %d}\n", i, j);
    VectorPair independent_pair = ComputeIndependentPair(i, j);
    vector add_part_ij = independent_pair.first;
    vector add_part_ji = independent_pair.second;
    // lock for i'th body
    pthread_mutex_lock(&(mutexes[i]));
    accelerations[i] = addVectors(accelerations[i], independent_pair.first);
    pthread_mutex_unlock(&(mutexes[i]));
    // lock for j'th body
    pthread_mutex_lock(&(mutexes[j]));
    accelerations[j] = addVectors(accelerations[j], independent_pair.second);
    pthread_mutex_unlock(&(mutexes[j]));
}

void resetAccelerations()
{
    int i;
    for (i = 0; i < bodies; i++)
    {
        accelerations[i].x = 0.f;
        accelerations[i].y = 0.f;
    }
}

void computeAccelerationsInParallel()
{
    int i;
    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS,NULL);
    for (i = 0; i < threads_num; i++) {
        int *arg = malloc(sizeof(*arg));
        *arg = i;
        if (pthread_create(&threads_pool[i], NULL, &runPairThread, (void*)arg) != 0) {
            perror("Failed to create the thread");
        }
    }
    for (i = 0; i < threads_num; i++) {
        if (pthread_join(threads_pool[i], NULL) != 0) {
            perror("Failed to join the thread");
        }
    }
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

void computeVelocitiy(int i)
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
        computeVelocitiy(i);
}

void computePositions()
{
    for(int i = 0; i < bodies; ++i)
        computePosition(i);
}

void simulate()
{
    resetAccelerations();
    computeAccelerationsInParallel();
    // computeAccelerations();
    computePositions();  // для каждого тела
    computeVelocities(); // для каждого тела
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
        simulate();
        for (j = 0; j < bodies; j++)
            printf("Body %d : %lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", j + 1, positions[j].x, positions[j].y, velocities[j].x, velocities[j].y, accelerations[j].x, accelerations[j].y);
    }
    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The elapsed time is %f seconds", time_spent);

    destructSystem();

    return 0;
}