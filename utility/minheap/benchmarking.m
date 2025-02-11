n_dim = 5;
n_nodes = 10;


% Time creation
timeit(@() create_list1(n_dim, n_nodes))
% timeit(@() create_list2(n_dim, n_nodes))
timeit(@() create_heap(n_dim, n_nodes))


% Create identical data
heap = minheap(n_nodes, @(a,b) a.LBFBOX - b.LBFBOX, boxdata);
box = boxdata(rand(n_dim, 2), randi(10000), rand(n_dim, 1), 0, 0, 0);
head1 = dlnode(box.BOX, box.LBFBOX, box.LBFPT, box.N2FLAG, box.N0FLAG, box.FLAG);
% head2 = dlnode2(box);
heap.add(box);
for i=1:n_nodes
    box = boxdata(rand(n_dim, 2), randi(10000), rand(n_dim, 1), 0, 0, 0);
    dlnode(box.BOX, box.LBFBOX, box.LBFPT, box.N2FLAG, box.N0FLAG, box.FLAG).insertAfter(head1);
%     dlnode2(box).insertAfter(head2);
    heap.add(box);
end


% Time extraction
timeit(@() extract1(head1))
% timeit(@() extract2(head2))
timeit(@() heap.extract())


function head = create_list1(n_dim, n_nodes)
    rand_dlnode = @() dlnode(rand(n_dim, 2), randi(10000), rand(n_dim, 1), 0, 0, 0);
    
    head = rand_dlnode();
    for i=1:n_nodes-1
        rand_dlnode().insertAfter(head);
    end
end

function head = create_list2(n_dim, n_nodes)
    rand_box = @() boxdata(rand(n_dim, 2), randi(10000), rand(n_dim, 1), 0, 0, 0);
    
    head = dlnode2(rand_box());
    for i=1:n_nodes-1
        dlnode2(rand_box()).insertAfter(head);
    end
end

function heap = create_heap(n_dim, n_nodes)
    rand_box = @() boxdata(rand(n_dim, 2), randi(10000), rand(n_dim, 1), 0, 0, 0);
    
    heap = minheap(n_nodes, @(a,b) a.LBFBOX - b.LBFBOX, boxdata);
    for i=1:n_nodes
        heap.add(rand_box());
    end
end


function node = extract1(head)
    node = head;

    currNode = node.Next;
    while ~isempty(currNode)
        if currNode.LBFBOX < node.LBFBOX
            node = currNode;
        end
        currNode = currNode.Next;
    end

    node.removeNode();
end


function node = extract2(head)
    node = head;

    currNode = node.Next;
    while ~isempty(currNode)
        if currNode.Data.LBFBOX > node.Data.LBFBOX
            node = currNode;
        end
        currNode = currNode.Next;
    end

    node.removeNode();
end
