classdef dlnode < handle
   % dlnode A class to represent a doubly-linked node.
   % Link multiple dlnode objects together to create linked lists.
   properties
      BOX;    % array of unsigned integers uint8
      LBFBOX; % floating point value to save lb f
      LBFPT;  % floating point array to save the point x, such that f(x)=lb f(X)
      N2FLAG; % no. of 2 flags in the box
      N0FLAG; % no. of 0 flags in the box
      FLAG;   % only used for ag76,ag78 to keep track of whether the box is marked for deletion for a particular tmax value or not
      sumFLAG; % only used for ag76,ag78 to store sum of the entries in FLAG
   end
   properties(SetAccess = private)
      Next = dlnode.empty
      Prev = dlnode.empty
   end
   
   methods
      function node = dlnode(BOX,LBFBOX,LBFPT,N2FLAG,N0FLAG,FLAG,sumFLAG)
         % Construct a dlnode object
         if nargin > 0
            node.BOX = BOX;
            node.LBFBOX=LBFBOX;
            node.LBFPT=LBFPT;
            node.N2FLAG=N2FLAG;
            node.N0FLAG=N0FLAG;
            node.FLAG=FLAG;
            node.sumFLAG=sumFLAG;
         end
      end
      
      function insertAfter(newNode, nodeBefore)
         % Insert newNode after nodeBefore.
         removeNode(newNode);
         newNode.Next = nodeBefore.Next;
         newNode.Prev = nodeBefore;
         if ~isempty(nodeBefore.Next)
            nodeBefore.Next.Prev = newNode;
         end
         nodeBefore.Next = newNode;
      end
      
      function insertBefore(newNode, nodeAfter)
         % Insert newNode before nodeAfter.
         removeNode(newNode);
         newNode.Next = nodeAfter;
         newNode.Prev = nodeAfter.Prev;
         if ~isempty(nodeAfter.Prev)
            nodeAfter.Prev.Next = newNode;
         end
         nodeAfter.Prev = newNode;
      end
      
      function removeNode(node)
         % Remove a node from a linked list.
         if ~isscalar(node)
            error('Input must be scalar')
         end
         prevNode = node.Prev;
         nextNode = node.Next;
         if ~isempty(prevNode)
            prevNode.Next = nextNode;
         end
         if ~isempty(nextNode)
            nextNode.Prev = prevNode;
         end
         node.Next = dlnode.empty;
         node.Prev = dlnode.empty;
      end
      
      function clearList(node)
         % Clear the list before
         % clearing list variable
         prev = node.Prev;
         next = node.Next;
         removeNode(node)
         while ~isempty(next)
            node = next;
            next = node.Next;
            removeNode(node);
         end
         while ~isempty(prev)
            node = prev;
            prev = node.Prev;
            removeNode(node)
         end
      end
   end
   
   methods (Access = private)
      function delete(node)
         clearList(node)
      end
   end
end