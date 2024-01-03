// destructor of iteration_hist class to make it a doubly-linked list to a singly-linked list
template <class myT>
iteration_hist<myT>::~iteration_hist<myT>()
{
    delete next;
}