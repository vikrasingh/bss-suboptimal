classdef minheap_till25Jan23 < handle
    % used till 25Jan23
    % Minimum binary heap, implemented using a dynamic array

    properties (Access = private)
        Capacity uint64 = 4  % current length of the dynamic array
        Array = zeros(1, 4)  % dynamic array
        Compare function_handle = @(a,b) a-b  % anonymous function used to compare the priorities of different elements
        Default = 0  % garbage value used only to store the type of the array
    end

    properties (SetAccess = private)
        Size uint64 = 0  % The current number of elements stored in the heap
    end
    
    methods
        function heap = minheap(capacity, compare, dtype)
            % Construct an empty binary heap with given capacity, compare function, and data type
            % Data type defaults to double, as represented by garbage value 0

            if nargin > 0
                heap.Capacity = capacity;
                if nargin > 1
                    heap.Compare = compare;
                end
                if nargin > 2
                    heap.Default = dtype;
                    a(1, capacity) = dtype;
                    heap.Array = a;
                else
                    heap.Array = zeros(1, capacity);
                end
            end
        end

        function add(heap, value)
            % Add a new value to the heap

            if heap.Size == heap.Capacity
                % Increase capacity, if no space left
                heap.doubleCapacity();
            end

            % Place new value at end of tree, and swap with lower-priority parents
            heap.Size = heap.Size + 1;
            heap.Array(heap.Size) = value;
            heap.siftUp(heap.Size);
        end

        function data = peek(heap)
            % Get highest-priority value (root), without removing it
            data = heap.Array(1);
        end

        function data = addFLAG(heap)
            % VS 30Dec23: add the FLAG variable for all the boxes in the heap
           data=[];
           if heap.Size>0
              data=heap.Array(1).FLAG; 
              if heap.Size>1
                  for i=2:heap.Size
                     data=data+heap.Array(i).FLAG;   
                  end
              end
           end
        end

        function data = extract(heap)
            % Get and remove highest-priority value (root)

            data = heap.Array(1);

            % Take last element, replace at root of heap, and swap with higher-priority children
            heap.Array(1) = heap.Array(heap.Size);
            heap.Size = heap.Size - 1;
            heap.siftDown(1);

            if heap.Capacity > 4 && heap.Size*4 <= heap.Capacity
                % Decrease capacity to free up memory
                heap.halfCapacity();
            end
        end
    end

    methods (Access = private)
        function heap = doubleCapacity(heap)
            % Double size of dynamic array
            heap.Capacity = 2*heap.Capacity;
            heap.Array(1, heap.Capacity) = heap.Default;  % This placing the garbage value at the desired end of the dynamic array
        end

        function heap = halfCapacity(heap)
            % Half size of dynamic array
            heap.Capacity = idivide(heap.Capacity, uint64(2));
            heap.Array = heap.Array(1:heap.Capacity);
        end

        function heap = siftUp(heap, index)
            % Repeatedly swap a node with lower-priority parent

            value = heap.Array(index);
            currI = index;
            while currI ~= 1 && heap.Compare(value, heap.Array(idivide(currI, uint64(2)))) < 0
                % While not at root and while parent is lower priority, make the swap
                heap.Array(currI) = heap.Array(idivide(currI, uint64(2)));
                currI = idivide(currI, uint64(2));
            end
            heap.Array(currI) = value;
        end

        function heap = siftDown(heap, index)
            % Repeatedly swap a node with higher-priority child

            value = heap.Array(index);
            currI = index;
            while true
                if 2*currI > heap.Size
                    % If leaf node, stop
                    break
                elseif heap.Compare(value, heap.Array(2*currI)) > 0
                    % If left child is higher-priority...
                    if 2*currI+1 <= heap.Size ...
                            && heap.Compare(heap.Array(2*currI), heap.Array(2*currI+1)) > 0
                        % ...and right child is even higher, swap with right child
                        heap.Array(currI) = heap.Array(2*currI+1);
                        currI = 2*currI+1;
                    else
                        % otherwise, swap with left child
                        heap.Array(currI) = heap.Array(2*currI);
                        currI = 2*currI;
                    end
                elseif 2*currI+1 <= heap.Size ...
                        && heap.Compare(value, heap.Array(2*currI+1)) > 0
                    % If right child is higher-priority, swap with right child
                    heap.Array(currI) = heap.Array(2*currI+1);
                    currI = 2*currI+1;
                else
                    % If no children are higher-priority, stop
                    break
                end
            end
            heap.Array(currI) = value;
        end
    end
end
