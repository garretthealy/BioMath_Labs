%% Lab 4 BIOS 26211: Line 'em up! Sequence alignment using dynamic programming 
%
% Name: Garrett Healy
%% Part 1. Global alignment using the Needleman-Wunsch algorithm
%
% Follow the pseudocode in the link here
% https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm to
% implement the N-W algorithm. 
% There are two separate parts of the algorithm: 1) filling the alignment
% matrix with scores and traceback directions, and 2) tracing back the
% alignment starting from the last two characters and printing out the
% result. You will implement those in two separate functions and test them
% before using them to align real protein sequences. 
%
% 1.1 N-W alignment scoring matrix
% First, write a function to caclulate the alignment scores and
% traceback matrix. It should receive three inputs: the two sequences (as
% strings/vectors of characters) and the gap penalty. 
% Instead of a scoring matrix, use binary scoring: +1 for a match
% and 0 for a mismatch. Within the function, I recommend using two separate
% matrices: 1) scoring matrix and 2) the traceback matrix, in which every cell
% contains an indication of which of the three options (insert, delete,
% match) resulted in the score for that cell. Remember that the matrices
% need an extra row and column to accommodate gaps before the first letter
% in either sequence. The function should output the traceback matrix and
% the optimal alignment score. 
%
% To test the function, calculate the alignment for two sequences given below. 
% Set the gap penalty of zero (gaps and mismatches are free, maximize
% matches!) and print out the resulting traceback matrix to check that it
% looks reasonable. The correct optimal alignment score is 5.
% Then set the gap penalty to 1, report the effect on the alignment matrix
% and the score, and print out the traceback matrix to check that it looks
% reasonable. The correct optimal alignment score is 0.
%
% seq1 = 'PAWHEAE'
% seq2 = 'HEAGAWGHEE'
%
% 1.2 Traceback of the alignment
% Implement the second half of the N-W algorithm. The function should take
% in three inputs: the traceback matrix (with characters or numbers
% representing the optimal choice made in every cell) and the two
% sequences. It should return the alignment as a three-row matrix, with the
% first row containing the first sequence with gaps indicated by dashes '-', 
% the second row containind '|' for every identical match and spaces
% everywhere else, and the second sequences in the third row.
%
% Test the function using the same two pairs of test sequences, again with
% gap penalty of 0 and then 1. The correct alignment with gap penalty 0 is:
% ---PAW-HEAE
%     || || |
% HEAGAWGHE-E
% The correct alignment with gap penalty 1 is:
% ---PAWHEAE
%     ||   |
% HEAGAWGHEE
%

clear all;
close all;

seq1 = 'PAWHEAE';
seq2 = 'HEAGAWGHEE';
gappen = 0;
mmpen = 0;

[scorematTest,tracematTest,alignvalTest,alignstrTest] = seqalign(seq1,seq2,gappen,mmpen);

gappen = 1;

[scorematTest2,tracematTest2,alignvalTest2,alignstrTest2] = seqalign(seq1,seq2,gappen,mmpen);
%
% 1.3 Application to Hemoglobin sequences
% Now you can use your functions to run alignments on actual sequences. In
% the file hemoglobins.mat you will find four character variables
% containing the following protein sequences: human hemoglobin alpha (humanHA),
% human hemoglobin beta (humanHB), horse hemoglobin beta (horseHB), and
% globin from the European river lamprey (lampreyH). Using your functions,
% align all six possible pairs of sequences using gap penalties of 1 and 2
% and report the number of idenitical letters in each one. Go to BLAST
% http://blast.ncbi.nlm.nih.gov/Blast.cgi select protein BLAST, click align
% 2 or more sequences, and align all six pairs of the four hemoglobin
% sequences (use the following accession numbers: P69905.2, P68871.2,
% P02062.1, P02207.2). Report and compare the identity fractions with the
% ones obtained by your own code. Which gap penalty agrees best with BLAST?
%

load('hemoglobins.mat');

gappen = 1;

seq1 = horseHB; 
seq2 = humanHA; 
[scoremat1,tracemat1,alignval1,alignstr1] = seqalign(seq1,seq2,gappen,mmpen);
count1 = length(find(alignstr1(2,:)));

seq1 = horseHB; 
seq2 = humanHB; 
[scoremat2,tracemat2,alignval2,alignstr2] = seqalign(seq1,seq2,gappen,mmpen);
count2 = length(find(alignstr2(2,:)));

seq1 = horseHB; 
seq2 = lampreyH; 
[scoremat3,tracemat3,alignval3,alignstr3] = seqalign(seq1,seq2,gappen,mmpen);
count3 = length(find(alignstr3(2,:)));

seq1 = humanHB; 
seq2 = humanHA; 
[scoremat4,tracemat4,alignval4,alignstr4] = seqalign(seq1,seq2,gappen,mmpen);
count4 = length(find(alignstr4(2,:)));

seq1 = lampreyH; 
seq2 = humanHA; 
[scoremat5,tracemat5,alignval5,alignstr5] = seqalign(seq1,seq2,gappen,mmpen);
count5 = length(find(alignstr5(2,:)));

seq1 = lampreyH; 
seq2 = humanHB; 
[scoremat6,tracemat6,alignval6,alignstr6] = seqalign(seq1,seq2,gappen,mmpen);
count6 = length(find(alignstr6(2,:)));

counts = zeros(2,6);
counts(1,:) = [count1 count2 count3 count4 count5 count6];

gappen = 2;

seq1 = horseHB; 
seq2 = humanHA; 
[scoremat1,tracemat1,alignval1,alignstr1] = seqalign(seq1,seq2,gappen,mmpen);
counts(2,1) = length(find(alignstr1(2,:)));

seq1 = horseHB; 
seq2 = humanHB; 
[scoremat2,tracemat2,alignval2,alignstr2] = seqalign(seq1,seq2,gappen,mmpen);
counts(2,2) = length(find(alignstr2(2,:)));

seq1 = horseHB; 
seq2 = lampreyH; 
[scoremat3,tracemat3,alignval3,alignstr3] = seqalign(seq1,seq2,gappen,mmpen);
counts(2,3) = length(find(alignstr3(2,:)));

seq1 = humanHB; 
seq2 = humanHA; 
[scoremat4,tracemat4,alignval4,alignstr4] = seqalign(seq1,seq2,gappen,mmpen);
counts(2,4) = length(find(alignstr4(2,:)));

seq1 = lampreyH; 
seq2 = humanHA; 
[scoremat5,tracemat5,alignval5,alignstr5] = seqalign(seq1,seq2,gappen,mmpen);
counts(2,5) = length(find(alignstr5(2,:)));

seq1 = lampreyH; 
seq2 = humanHB; 
[scoremat6,tracemat6,alignval6,alignstr6] = seqalign(seq1,seq2,gappen,mmpen);
counts(2,6) = length(find(alignstr6(2,:)));

%% 
% 
% It appears that BLAST uses a gap penalty close to 1, since it's values
% are very similar to the values that were obtained with a gap penalty of 1
% , as opposed to those obtained with a penalty of 2. 
%

%% Part 2: Local alignment using the Smith-Waterman algorithm
%
% 2.1 S-W scoring matrix
% Write a function to calculate the scoring matrix following the
% Smith-Waterman algobrithm, which can be done by modifying the N-W function;
% use tne score of 1 for a match and 0 for mismatch. Input: two
% sequences and the gap penalty; output: the traceback matrix
% and the scoring matrix (since the traceback starts at the maximum value
% in the scorin matrix.)
%
% To test the function, align the short sequence 'PNLHGLFGRKTG' with the
% Cytochrome C sequence found in variable CytC in the data file
% sequences.mat; use gap penalty of 2. The sequence is a motif within
% Cytochrome Cs so the maximum score in the scoring matrix should be 12. 
%
% 2.2. S-W alignment traceback
% Write a function to output the alignment based the Smith-Waterman
% traceback matrix (or the scoring matrix); this requires only a minor
% modification of the N-W alignment function from part 1. Input:
% scoring matrix, traceback matrix, two sequences; output: local alignment
% in the same format as in part 1. Test your function by aligning the sequence
% 'PNLHGLFGRKTG' with the Cytochrome C sequence found in variable CytC in
% the data file sequences.mat; use gap penalty of 2. The alignment should
% show all 12 letters of the short sequence perfectly aligned with the
% corresponding subsequence in CytC.
%

clear all;
close all;

load('sequences.mat');

seq1 = 'PNLHGLFGRKTG';
seq2 = CytC;
gappen = 2;
mmpen = 0;

[scoremat,tracemat,alignval,alignstr] = swseqalign(seq1,seq2,gappen,mmpen);
alignvalTest = alignval;
alignstrTest = alignstr;

% 2.3 Finding similarity in distant sequences
% In this exercise you will use local alignment to compare distantly
% related protein sequences from different organisms saved in the file
% sequences.mat. Align the sequences RBP (retinol-binding protein 4) with
% BBP (bilin-binding protein) using gap penalty of 2. Compare the results
% with the results returned by BLAST (use accession numbers NP_006735.2 and
% P09464.2). Align the sequences BBP with APO (Apolipoprotein)  using the
% gap penalty of 2. Compare the results with the results returned by BLAST
% (use accession numbers P05090.1 and P09464.2). Try changing the gap
% penalty to 1 and 3 and report which one results in number of idenitical
% letters closest to BLAST. 


[scoremat,tracemat,alignval,alignstr] = swseqalign(RBP,BBP,gappen,mmpen);
alignvalRBP_BBP2 = alignval;
countRBP_BBP2 = length(find(alignstr(2,:)));

[scoremat,tracemat,alignval,alignstr] = swseqalign(APO,BBP,gappen,mmpen);
alignvalAPO_BBP2 = alignval;
countAPO_BBP2 = length(find(alignstr(2,:)));
gappen = 1;
[scoremat,tracemat,alignval,alignstr] = swseqalign(RBP,BBP,gappen,mmpen);
alignvalRBP_BBP1 = alignval;
alignstrRBP_BBP1 = alignstr;
countRBP_BBP1 = length(find(alignstr(2,:)));

[scoremat,tracemat,alignval,alignstr] = swseqalign(APO,BBP,gappen,mmpen);
alignvalAPO_BBP1 = alignval;
countAPO_BBP1 = length(find(alignstr(2,:)));
gappen = 3;
[scoremat,tracemat,alignval,alignstr] = swseqalign(RBP,BBP,gappen,mmpen);
alignvalRBP_BBP3 = alignval;
countRBP_BBP3 = length(find(alignstr(2,:)));
[scoremat,tracemat,alignval,alignstr] = swseqalign(APO,BBP,gappen,mmpen);
alignvalAPO_BBP3 = alignval;
countAPO_BBP3 = length(find(alignstr(2,:)));


%%
% By looking at BLAST, it can be seen that it is more closely related to
% the gap penalty of 3 when using the Smith-Waterman algorithm.
