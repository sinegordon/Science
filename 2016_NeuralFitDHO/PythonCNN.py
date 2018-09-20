
# coding: utf-8

# In[1]:

import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn import model_selection, metrics, datasets
from neupy import algorithms, layers, environment


environment.reproducible()
environment.speedup()


# In[32]:

nmax = 4

def dhofun(x, nl, *p):
    res = 0.0    
    for i in xrange(nl):
        A = p[3*i]
        G = p[3*i+1]
        E = p[3*i+2]
        X2 = x**2
        E2 = E*E
        res += X2*E2*G*A/((X2 - E2)**2 + X2*G*G)
    maxl = np.max(res)
    return res/maxl

def get_data(count):
    
    x = np.linspace(1, 10, 100)
    x_train = np.zeros((count, 1, 1, 100))
    x_test = np.zeros((100, 1, 1, 100))
    y_train = np.zeros((count, nmax))
    y_test = np.zeros((100, nmax))
    p = np.zeros(3*nmax)
    for i in range(count):
        nl = np.random.randint(1, nmax+1)
        for k in range(nl):
            p[3*k] = np.random.uniform(0, 1)
            p[3*k+1] = np.random.uniform(0, 5)
            p[3*k+2] = np.random.uniform(1, 10)
        x_train[i, 0, 0, :] = dhofun(x, nl, *p)
        y_train[i, nl-1] = 1

    for i in range(100):
        nl = np.random.randint(1, nmax+1)
        for k in range(nl):
            p[3*k] = np.random.uniform(0, 1)
            p[3*k+1] = np.random.uniform(0, 5)
            p[3*k+2] = np.random.uniform(1, 10)
        x_test[i, 0, 0, :] = dhofun(x, nl, *p)
        y_test[i, nl-1] = 1
               
    return x_train, x_test, y_train, y_test


network = algorithms.Adadelta(
    [
        layers.Input((1, 1, 100)),

        layers.Convolution((16, 1, 7)) > layers.BatchNorm() > layers.Relu(),
        layers.MaxPooling((1, 3)),
        layers.Convolution((8, 1, 5)) > layers.BatchNorm() > layers.Relu(),
        layers.MaxPooling((1, 2)),
        layers.Convolution((4, 1, 3)) > layers.BatchNorm() > layers.Relu(),
        layers.MaxPooling((1, 2)),

        layers.Reshape(),
        layers.Linear(1024) > layers.BatchNorm() > layers.Relu(),
        layers.Softmax(nmax),
    ],

    # Using categorical cross-entropy as a loss function
    error='categorical_crossentropy',

    # Min-batch size
    batch_size=128,

    # Learning rate. We can allow high values
    # since we are using Batch Normalization
    step=1.0,

    # Shows information about algorithm and
    # training progress in terminal
    verbose=True,

    # Randomly shuffles training dataset before every epoch
    shuffle_data=True,

    # Step decay algorithm minimizes learning step
    # monotonically after each iteration.
    addons=[algorithms.StepDecay],
    # Parameter controls step redution frequency. The higher
    # the value the slower step parameter decreases.
    reduction_freq=8,
)


# In[33]:

# Shows networks architecture in terminal's output
network.architecture()

print("Generate data...")
x_train, x_test, y_train, y_test = get_data(100000)
print("Done generate data.")


# In[34]:

# Train for only two epochs
network.train(x_train, y_train, x_test, y_test, epochs=4)


# In[37]:

y_predicted = network.predict(x_test).argmax(axis=1)
y_test_labels = np.asarray(y_test.argmax(axis=1)).reshape(len(y_test))

print(metrics.classification_report(y_test_labels, y_predicted))
score = metrics.accuracy_score(y_test_labels, y_predicted)
print("Validation accuracy: {:.2%}".format(score))


# In[ ]:



