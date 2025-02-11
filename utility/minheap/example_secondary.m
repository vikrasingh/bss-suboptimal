n_dim = 20;
n_nodes = 100000;


% Create data.
minimize_LBFBOX = @(a,b) a.LBFBOX - b.LBFBOX;
% heap = minkheap(3, n_nodes, minimize_LBFBOX, boxdata);
heap = minheap( n_nodes, minimize_LBFBOX, boxdata);
for i=1:n_nodes
    heap.add(boxdata(rand(n_dim, 2), 10000*rand(1), rand(n_dim, 1), randi(n_dim), 0, 0,[]));
end


% Extract minimum for primary criterion.
box = heap.extract();
box.LBFBOX
box.N2FLAG


% Extract estimated minimum for secondary criterion.
% Sample is first 10 (approximately the 10 best) and 40 random elements.
% Note that we wish to maximize N2FLAG, so the compare function should give
% a smaller value when a.N2FLAG > b.N2FLAG.
maximize_N2FLAG = @(a,b) b.N2FLAG - a.N2FLAG;
box = heap.extract_secondary(10, 40, maximize_N2FLAG);
box.LBFBOX
box.N2FLAG


% extract_secondary works in the same way for minheap
