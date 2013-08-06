## NAME

    BoundaryFinder - Exon/Intron Boundary Finder Module

## SYNOPSIS

      use BoundaryFinder;

## DESCRIPTION

    The BoundaryFinder uses consensus sequences to find and score
    boundaries in a sequence. This was developed for finding exon
    boundaries, but it can be used with any consensus sequence.

    The boundary finder works best when you finding boundaries in 
    small sequences. This is because the number of matches is often
    large, especially if there are no "100%" positions. (e.g. 100% 
    of nucleotides have to be an "A" at position x, or it does not 
    match)

## OVERVIEW

    The general procedure of using the library to find exon boundaries
    is as follows:

### 1. Use the BlastManager to blast your query sequence agaist the
          subject database.

          This will give a BlastResults oject, or a set of blast hits.

### 2. Determine what BlastResults or hits you want use.

          Sometimes blast returns hits that are not relevant to your
          query. These could be weak matches that blast returned that
          do not make sense with your query. You can trim these out 
          by changing the blast parameters (such as a minimum "e" value)
          or iterating through the BlastResults and removing the 
          unwanted hits.

### 3. Use BoundaryFinder to find the boundaries of the blast hits

          Each query is a pair of sequences representing the 5` and 3`
          ends of a blast hit.

          This can be automatically done by passing a BoundaryResults 
          object, or done manually. It will create a query with a 50 
          extra nucleotides from each end of the edges of the blast hit. 
          You can easily specify a buffer size other than 50.

          ***************************************************************
          E.g. 

          This hastily done ascii example shows how the BoundaryFinder will
          create queries automatically from a hit. Given the blast result below,
          a sequence is created by using 50 nucleotides upstream (5') to the end
          of the blast hit (3' end). This is also done for the 3' end, with 
	  50 NT after the end of the hit.

          Query   ATGATGATGATGATG
                  | |||||||||||    ... and so on
          Sub     ACGATGATGATGACC
                 
          |< -50  ^  ------------- To the end of the blast hit --------------->|

          ^ = Points to the 5' edge of the blast result.

          ***************************************************************

          This can be done manually too, and it takes only a little more 
          work. You will just have to figure out the buffers for each hit,
          and add them with the BoundaryFinder::addQuery method.

### 4. Run the BoundaryFinder to get a BoundaryResults object.

          The BoundaryResults object contains all the non-zero scores, and
          their positions in the blast database for easy retrieval. It is 
          simply a wrapper of two hashes that contain the position as keys,
          and score as the value.

## USAGE

     Coming Soon...

## AUTHOR

        Nathan Elmore
        nate@elmoren.com
        http://www.elmoren.com

## COPYRIGHT

	

## SEE ALSO
   
	Blast+ 2.2.28

