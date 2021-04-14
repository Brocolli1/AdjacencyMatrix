#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>
#include <math.h>
#include <string.h>

using namespace std;


typedef struct {
    PyObject_HEAD
    long long vertices;
    long long matrix[64];
} AdjacencyMatrix;

static PyObject *AdjacencyMatrix__new__(PyTypeObject *type, PyObject *args) {
    return type->tp_alloc(type, 0);
}

static void AdjacencyMatrix__del__(AdjacencyMatrix *self) {
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static int AdjacencyMatrix__init__(AdjacencyMatrix *self, PyObject *args) {
    const char *graph;
    char a; 
    int length; 
    for(int i = 0; i < 64; i++)
    {
        self->matrix[i] = 0; 
    }
    if(PyArg_ParseTuple(args, "s#", &graph, &length))
    {

        self->vertices = 0;
        for(int i = 0; i < graph[0]-63; i++)
        {
            _int64 mask = (_int64)1 << i; 
            self->vertices = self->vertices | mask; 
        }
        int position = 0, line = 1; 
        for(int i = 1; i < length; i++){
            a = graph[i] - 63;
            a << 2; 
            do
            {
                if(a & (1 << 7))
                    self->matrix[line] = self->matrix[line] | (1 << (63 - position));
                position++;
                if(position == line)
                { 
                    position = 0;
                    line++;
                }
            } while (a = a << 1);
        }
    }
    else {
        self->vertices = 0; 
    }
    return 0;
}

static PyObject *number_of_vertices(AdjacencyMatrix *self) {
    _int64 vertices = self->vertices;
    int result = 0;
    do 
    {
        result += vertices % 2;
    } while( vertices = vertices >> 1 );

    return PyLong_FromLong(result); 
}

static PyObject *vertices(AdjacencyMatrix *self) {
    _int64 vertices = self->vertices;
    PyObject*set = PySet_New(nullptr);
    int i = 0; 
    do 
    {
        if(vertices % 2)
        PySet_Add(set, PyLong_FromLong(i));
        i++; 
    } while( vertices = vertices >> 1 );

    return set;
}

static PyObject *vertex_degree( AdjacencyMatrix *self, PyObject*args)
{
    int vertex; 
    PyArg_ParseTuple(args, "i", &vertex);
    _int64 neighbours = self->matrix[vertex];
    int result = 0;
    do 
    {
        result += neighbours % 2;
    } while( neighbours = neighbours >> 1 );

    return PyLong_FromLong(result); 
} 

static PyObject *vertex_neighbors( AdjacencyMatrix *self, PyObject*args)
{
    int vertex; 
    PyArg_ParseTuple(args, "i", &vertex);
    _int64 neighbours = self->matrix[vertex];
    PyObject*set = PySet_New(nullptr);
    int i = 1; 
    do 
    {
        if(neighbours % 2)
        PySet_Add(set, PyLong_FromLong(i));
        i++; 
    } while( neighbours = neighbours >> 1 );

    return set;
} 

static PyObject *add_vertex( AdjacencyMatrix *self, PyObject*args)
{
    int v; 
    PyArg_ParseTuple(args, "i", &v);
    if(v<64)
    {
        
        _int64 mask = (_int64) 1 << v;  
        self->vertices = self->vertices | mask;  
    }
    Py_RETURN_NONE;
}

static PyObject *delete_vertex( AdjacencyMatrix *self, PyObject*args)
{
    int v; 
    PyArg_ParseTuple(args, "i", &v);
    if(v<64)
    {
        self->matrix[v] = 0;
        _int64 mask = (~0) ^ (1 << v); 
        self->vertices = self->vertices & mask;
        for(int i = 0; i < 64; i++)
        {
            self->matrix[v] = self->matrix[v] & mask;
        }
        
    }
    Py_RETURN_NONE;
}

static PyObject *number_of_edges( AdjacencyMatrix *self)
    {
        int sum = 0;
        for(int i = 1; i < 64; i++)
        {
            int j =0;
            _int64 neighbours = self->matrix[i];
            int result = 0;
            for(int j = 0; j < i; j++ );
            {
                sum += neighbours % 2;
                neighbours = neighbours >> 1;
                j++;
            }
        }
        return PyLong_FromLong(sum);
    }

static PyObject *edges(AdjacencyMatrix *self)
    {
        _int64 vertices = self->vertices;
        PyObject*set = PySet_New(nullptr);
        int i = 1; 
        
        int position = 0, line = 1; 
        do
            {
                
                if(self->matrix[line] & (1 << (63 - position)))
                {
                    PyObject*edge;
                    edge = PyTuple_New(2);
                    PyTuple_SetItem(edge, 0, PyLong_FromLong(position));
                    PyTuple_SetItem(edge, 1, PyLong_FromLong(line));
                    PySet_Add(set, edge);
                }
                position++;
                if(position == line)
                { 
                    position = 0;
                    line++;
                }
            } while (line < 64);    
        return set;
    }

 static PyObject *is_edge(AdjacencyMatrix *self, PyObject*args)
    {
    int v, u; 
    PyArg_ParseTuple(args, "ii", &v, &u);
        _int64 mask = (~0) ^ (1 << v); 
        return PyBool_FromLong(self->matrix[u] & mask);
    }

static PyObject *add_edge( AdjacencyMatrix *self,  PyObject*args)
{
    int v, u; 
    PyArg_ParseTuple(args, "ii", &v, &u);
    _int64 mask1 = 1 << (63-v); 
    _int64 mask2 = 1 << (63-u);
    self->matrix[u] = self->matrix[u] | mask1;
    self->matrix[v] = self->matrix[v] | mask2;  
     Py_RETURN_NONE;
}

static PyObject *delete_edge( AdjacencyMatrix *self,  PyObject*args)
{
    int v, u; 
    PyArg_ParseTuple(args, "ii", &v, &u);
    _int64 mask1 = (~0) ^ (1 << (63-v)); 
    _int64 mask2 = (~0) ^ (1 << (63-u)); 
    self->matrix[u] = self->matrix[u] & mask1;
    self->matrix[v] = self->matrix[v] & mask2;
     Py_RETURN_NONE;
}

static PyObject *is_complete_bipartite( AdjacencyMatrix *self)
{
    for(int i = 1; i<64; i++)
    {
        _int64 mask = (~0) ^ (1 << i); 
        for(int j = 0; j<i; j++)
        {
            if(self->matrix[j] & mask)
            {
                if(self->matrix[i] & self->matrix[j]) 
                    return PyBool_FromLong(0);
            }
        }
    } 
    return PyBool_FromLong(1);
}

static PyObject *
Largest_richcompare(AdjacencyMatrix *self, AdjacencyMatrix *other, int op)
{
    PyObject *result = NULL;
    result = Py_True;
    if(self->vertices != other->vertices) 
    {
        result = Py_False;        
    }
    for(int i = 0; i < 64; i++)
    {
        if(self->matrix[i] != other->matrix[i])
        {
            result = Py_False;
            break;
        }
    }     

    switch (op) {
        case Py_EQ:
            Py_XINCREF(result);
            return result;
            break; 
        case Py_NE:
            if(result == Py_False)
                result = Py_True;
            else result = Py_False;
            Py_XINCREF(result);
            return result;
            break; 
    }
}

static PyMethodDef AdjacencyMatrixMethods[] = {
    { "number_of_vertices", (PyCFunction)number_of_vertices, METH_NOARGS, "Returns number of vertices."},
    { "vertices", (PyCFunction)vertices, METH_NOARGS, "Returns a set of all vertices."},
    { "vertex_degree", (PyCFunction)vertex_degree, METH_VARARGS, "Returns degree of a given vertex."},
    { "vertex_neighbors", (PyCFunction)vertex_neighbors, METH_VARARGS, "Returns degree of a given vertex."},
    { "delete_vertex", (PyCFunction)delete_vertex, METH_VARARGS, "Deletes a given vertex."},
    { "add_vertex", (PyCFunction)add_vertex, METH_VARARGS, "Adds a given vertex."},
    { "delete_vertex", (PyCFunction)delete_vertex, METH_VARARGS, "Deletes a given vertex."},
    { "number_of_edges", (PyCFunction)number_of_edges, METH_NOARGS, "Returns a number of all adges in the graph."},
    { "edges", (PyCFunction)edges, METH_NOARGS, "Returns a set all edges in a graph."},
    { "is_edge", (PyCFunction)is_edge, METH_VARARGS, "Checks if an edge exists between given vertexes."},
    { "add_edge", (PyCFunction)add_edge, METH_VARARGS, "Adds an adge beteen given vertexes."},
    { "delete_edge", (PyCFunction)delete_edge, METH_VARARGS, "Deletes an edge between given vertexes."},
    { "is_complete_bipartite", (PyCFunction)is_complete_bipartite, METH_NOARGS, "Checks if the graph is complete bipartite."},
    { NULL }
};


  



static PyTypeObject AdjacencyMatrixType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "AdjacencyMatrix.AdjacencyMatrix",
    sizeof(AdjacencyMatrix),
    0,
    (destructor)AdjacencyMatrix__del__,
    0,0,0,0,0,0,0,0,0,0,
    0,                          
    0,0,0,
    Py_TPFLAGS_DEFAULT,
    "Adjacency Matrix.",
    0,0,
    (richcmpfunc)&Largest_richcompare,
    0,0,0,
    AdjacencyMatrixMethods,
    0,0,0,0,0,0,0,
    (initproc)AdjacencyMatrix__init__,
    0,
    (newfunc)AdjacencyMatrix__new__
};

static PyModuleDef simple_graphs_module = {
    PyModuleDef_HEAD_INIT,
    "simple_graphs",
    "Adjacency Matrix.",
    -1,
    NULL, NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC PyInit_simple_graphs(void) {
    if (PyType_Ready( &AdjacencyMatrixType ) < 0) return NULL;

    PyObject* m = PyModule_Create( &simple_graphs_module );
    if (m == NULL) return NULL;
    
    Py_INCREF( &AdjacencyMatrixType );
    PyModule_AddObject( m, "AdjacencyMatrix",(PyObject *)&AdjacencyMatrixType);
    return m;
}


