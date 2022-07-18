#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('lsmagic', '')


# In[2]:


a = 5

get_ipython().run_line_magic('whos', '')


# In[3]:


get_ipython().show_usage()


# In[4]:


get_ipython().run_line_magic('quickref', '')


# In[5]:


get_ipython().run_line_magic('pinfo', 'round')


# In[6]:


a = 5
b = 5.0
c = float(5)
d = 'dee'
e = 'e'

type(a), type(b), type(c), type(d), type(e)


# In[7]:


empty_list = []

string_list = ['lions', 'tigers', 'bears', 'sharks', 'hamsters']

int_list = [0, 1, 2, 3, 4]

int_list2 = list(range(5,10))

list_from_variables = [a,b,c,d,e]

list_of_lists = [empty_list,
                 string_list,
                 list_from_variables,
                 int_list,
                 int_list2]


# In[8]:


empty_list


# In[9]:


string_list


# In[10]:


int_list


# In[11]:


int_list2


# In[12]:


list_from_variables


# In[13]:


list_of_lists


# In[14]:


string_list[0]


# In[15]:


string_list[1:4]


# In[16]:


int_list[::2]


# In[17]:


list_of_lists[1][4]


# In[18]:


int_list[2] = 222

int_list.append(5)

string_list.remove('lions')

list_from_variables.extend(int_list)


# In[19]:


int_list


# In[20]:


string_list


# In[21]:


list_from_variables


# In[22]:


joe_blow = (32, 'tall', 'likes hats')
joe_blow


# In[23]:


joe_blow[1]


# In[24]:


joe_blow.append('married')


# In[25]:


joe_blow[2] = 'not really a fan of hats'


# In[26]:


pets = ('elephant', 'cow', 'rock')

pet1, pet2, pet3 = pets

pet1


# In[27]:


tuple_with_one_item = pet1,

tuple_with_one_item


# In[28]:


dictionary0 = {'key1': 'value1', 'key2': 'value2', 'key3': 'value3'}
dictionary0


# In[29]:


dictionary1 = {1: 'value1', 2: 'value2', 3: 'value3'}
dictionary1


# In[30]:


list(dictionary1.keys())


# In[31]:


list(dictionary1.values())


# In[32]:


cylinder = {'mass': 50, 'base': 10, 'height': 100}
cylinder['mass']


# In[33]:


keys = ['mass01', 'inertia01', 'mass02', 'inertia02']
values = [10, 1, 50, 5]
dict(zip(keys, values))


# In[34]:


def abs_value(A):
    if A < 0:
        A = -A
    return A

abs_value(-100)


# In[35]:


abs_value(123)


# In[36]:


def long_div(dividend, divisor):
    quotient = dividend // divisor  # // : floor division
    remainder = dividend % divisor  # % : modulo
    return quotient, remainder


# In[37]:


a = 430
b = 25

quo, rem = long_div(a, b)

quo, rem


# In[38]:


msg = '{} divided {} is {} remainder {}'.format(a, b, quo, rem)
print(msg)


# In[39]:


import sys

print(sys.version)


# In[40]:


from sys import version

print(version)


# In[41]:


import sympy as sm
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


# In[42]:


sm.cos(12.0)


# In[43]:


np.cos(12.0)


# In[44]:


sp.cos(12.0)

