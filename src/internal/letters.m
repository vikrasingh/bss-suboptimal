% to get a string of excel columns A,B,...,Z,AA,...
function [alphabetCells]=letters

    letters=num2cell(char(65:90));   % to get 1x26 cell array containing letters from A-Z
    alphabetCells=cell(1,1378);  % so far  A,...,Z,AA,AB,...,AZ,BA,BB,...BZ,CA,CB,...CZ  1378 columns can be handled
    alphabetCells(1:26)=letters;
    for letterctr=2:27    % upper bnd 4 because only 4 set of 26 letters has been used
        for secondctr=1:26
            alphabetCells{ secondctr+26*(letterctr-1) }=append( letters{letterctr-1},letters{secondctr} ) ;
        end
    end

    for tripleletter=1:1   % 1: 1 because we are creating only AAA...AAZ, ABA...ABZ, . . . , AZA... AZZ   
       alphabetCells( 1+26+676*tripleletter : 26+676*(tripleletter+1) )=cellfun(@(c)strcat(letters{tripleletter}, c) , alphabetCells( 1+26+676*(tripleletter-1):26+676*tripleletter ) ,'UniformOutput',false);
    end

    
end