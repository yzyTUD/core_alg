#include <pthread.h>
#include <stdio.h>

#define PRODUCER_COUNT 1
#define CONSUMER_COUNT 1
#define THREAD_COUNT (PRODUCER_COUNT + CONSUMER_COUNT)

#define QUEUE_SIZE 64

// Let's go the easy way and keep a gap between head and tail when full.
static volatile int queue_head = 1;
static int queue[QUEUE_SIZE];
static volatile int queue_tail = 0;

static void *producer(void *data_ptr);
static void *consumer(void *data_ptr);

int main()
{
    pthread_t threads[THREAD_COUNT];
    int data[THREAD_COUNT];

    int thread_id;
    for(thread_id = 0; thread_id < PRODUCER_COUNT; thread_id++)
    {
        data[thread_id] = thread_id;
        pthread_create(&threads[thread_id], NULL, producer, &data[thread_id]);
    }

    for(; thread_id < THREAD_COUNT; thread_id++)
    {
        data[thread_id] = thread_id;
        pthread_create(&threads[thread_id], NULL, consumer, &data[thread_id]);
    }

    for(int i = 0; i < THREAD_COUNT; i++)
    {
        pthread_join(threads[i], NULL);
    }

    return 0;
}

static inline int advance(volatile int *idx)
{
    int old, new;
    do
    {
        old = *idx;
        new = (old + 1) % QUEUE_SIZE;
    } while(!__sync_bool_compare_and_swap(idx, old, new));
    return old;
}

static void *producer(void *data_ptr)
{
    int thread_id = *(int *)data_ptr;
    printf("[%d] spinning as producer\n", thread_id);

    for(int i = 0; i < 1000; i++)
    {
        while((queue_head + 1) % QUEUE_SIZE == queue_tail)
        {
            // spin like we've never spun before?
            sleep(0);
        }

        queue[queue_head] = i;
        advance(&queue_head);

        printf("[%d] produced %d\n", thread_id, i);
    }

    return NULL;
}

static void *consumer(void *data_ptr)
{
    int thread_id = *(int *)data_ptr;
    printf("[%d] consuming\n", thread_id);

    // instead of poison pill let's just consume exactly what is produced.
    for(int i = 0; i < 1000; i++)
    {
        while(queue_tail == queue_head)
        {
            // spin spin spin
            sleep(0);
        }

        int idx = advance(&queue_tail);
        printf("[%d] consumed %d\n", thread_id, queue[idx]);
    }

    printf("[%d] finished consuming\n", thread_id);
    return NULL;
}
