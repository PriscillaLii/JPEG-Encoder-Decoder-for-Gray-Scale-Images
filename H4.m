% Data Compression HW 4
% Writer: Sibei Li
% DEC/2017

% ----------     Read Data     ----------
[I,map]=imread('river.gif');
G=ind2gray(I,map);
Gdouble = double(G);

% ----------     Initialize Parameters     ----------
% Scalar
Qparam = 2.365;
% DPCM
DPCMparam = 1;

% ---------- Test Data(Not included in homework) ----------
%G = [62 55 55 54 49 48 47 55;
%62 57 54 52 48 47 48 53;
%61 60 52 49 48 47 49 54;
%63 61 60 60 63 65 68 65;
%67 67 70 74 79 85 91 92;
%82 95 101 106 114 115 112 117;
%96 111 115 119 128 128 130 127;
%109 121 127 133 139 141 140 133];

% Quantization matrix given in lecture node
Q = [16	11	10	16	24	40	51	61; 
12	12	14	19	26	58	60	55; 
14	13	16	24	40	57	69	56; 
14	17	22	29	51	87	80	62; 
18	22	37	56	68	109	103	77; 
24	35	55	64	81	104	113	92; 
49	64	78	87	103	121	120	101; 
72	92	95	98	112	100	103	99];
% ----------     Pretreatment: Nomalization, DCT, Divide by Q     ----------
[rows, cols] = size(Gdouble);
orgnSize = rows*cols*8;
newG = [];
for i = 1:rows/8
    newRow = [];
    for j = 1:cols/8
        temp = Gdouble( (i-1)*8+1 : i*8 , (j-1)*8+1 : j*8);
        % Mean-normalization
        temp = temp - 128;
        % Apply dct2 on every 8*8 block in G
        temp = dct2(temp);
        % Divide by Quantization Matrix
        temp = temp./(Qparam * Q);
        newRow = [newRow, temp];
    end
    newG = [newG; newRow];
end

Ground = round(newG);

% ----------     Entropy-coding of the DC coefficients     ----------
% *****   Extract DC coefficients   *****
DCs = [];
for i = 1:rows/8
    for j = 1:cols/8
        DCs = [DCs, Ground((i-1)*8+1 , (j-1)*8+1)];
    end
end
% *****   DPCM   *****
[DCrow, DCcol] = size(DCs);
DCdpcm = zeros(DCrow, DCcol);
DCdpcm(1) = DCs(1);
for i = 2:DCcol
    DCdpcm(i) = DCs(i) - DPCMparam * DCs(i-1);
end
% *****   Huffman code & Convert value into (hsm).   *****
% 1. find s
sizeDCdpcm = size(DCdpcm);
SofDC = zeros(1, sizeDCdpcm(2));
for i = 1 : sizeDCdpcm(2)
    if(DCdpcm(i) > 0 | DCdpcm(i) == 0)
        SofDC(i) = 1;
    end
end
% 2. use 2^(k-1) + t to represent the rest. Find k and t.
KofDC = zeros(1, sizeDCdpcm(2));
TofDC = zeros(1, sizeDCdpcm(2));
absDCdpcm = abs(DCdpcm);
for i = 1 : sizeDCdpcm(2)
    if (absDCdpcm(i) == 0)
        KofDC(i) = 0;
        TofDC(i) = 0;
    else
        KofDC(i) = floor(log2(absDCdpcm(i))) + 1;
        TofDC(i) = absDCdpcm(i) - pow2(KofDC(i) - 1);
    end
end
% 3. huffman code k
% (1). find the probability of each k.
CountOfK = zeros(1, 12);
for i = 1 : sizeDCdpcm(2)
    CountOfK(KofDC(i) + 1) = CountOfK(KofDC(i) + 1) + 1;
end
ProbofK = zeros(1, 12);
for i = 1 : 12
    ProbofK(i) = CountOfK(i) / sizeDCdpcm(2);
end
% (2). Huffman code
Hdict = huffmandict([0,1,2,3,4,5,6,7,8,9,10,11],ProbofK);
HofDC1 = cell(1, sizeDCdpcm(2));
for i = 1 : sizeDCdpcm(2)
    HofDC1{1, i} = huffmanenco(KofDC(i), Hdict);
end
% 4. Convert t into binary. 
% Use de2bi, because Mac OS cannot access Data Acquisition Toolbox, which contains function decimalToBinaryVector(). Notice it is inverse. 
MofDC = cell(1, sizeDCdpcm(2));
for i = 1 : sizeDCdpcm(2)
    if KofDC(i) == 1 | KofDC(i) == 0
        temp1 = [];
    else
        temp = de2bi(TofDC(i),KofDC(i)-1);
        [tempRow, tempCol] = size(temp);
        temp1 = [];
        for j = 0 : (tempCol - 1)
            temp1 = [temp1; temp(tempCol - j)];
        end
    end
    MofDC{1, i} = temp1;
end
% ----------     Now DC coefficients can be encoded.     ----------



% ----------     Entropy-coding of the AC coefficients     ----------
% Init some variables
ACs = [];
DofAC = [];
HofAC = [];
SofAC = [];
MofAC = [];
numOfBlock = 0;
% Extract each AC terms
% represent all non-zero items as (numOfZerosBefore, item.value)
for i = 1:rows/8
    for j = 1:cols/8
        numOfBlock = numOfBlock + 1;
        % block is current 8*8 block
        block = [Ground((i-1)*8+1 : (i-1)*8 + 8 , (j-1)*8+1: (j-1)*8+8)];
        % zigzag scan
        zzBlock = zigzag(block);
        % exclude the DC coefficients
        zzBlock = zzBlock(2: end);
        % represent each non-zero item in this block as (numOfZerosBefore, item.value)
        numOfZerosBefore = 0;
        ACnonZero = [];
        for i = 1 : 63
            if zzBlock(i) ~= 0
                temp = [numOfZerosBefore; zzBlock(i)];
                ACnonZero = [ACnonZero, temp];
                numOfZerosBefore = 0;
            else
                numOfZerosBefore = numOfZerosBefore + 1;
            end
            % Mark the end of the block
            if i == 63 & zzBlock(i) == 0
                temp = [0; 0];
                ACnonZero = [ACnonZero, temp];
            end
        end
        ACs = [ACs, ACnonZero];
    end
end
% d
DofAC = ACs(1,:);
% Now the first row is d (except every end of the block).
% *****   Huffman code & Convert item.value into (hsm).   *****
% 1. find s
ACsize = size(ACs);
SofAC = zeros(1,ACsize(2));
for i = 1 : ACsize(2)
    if ACs(1, i) == 0 & ACs(2, i) == 0
        % Mark the end of the block
        SofAC(i) = -1;
        continue;
    end
    if(ACs(2, i) > 0 | ACs(2, i) == 0)
        SofAC(i) = 1;
    end
end
% 2. use 2^(k-1) + t to represent the rest. Find k and t.
KofAC = zeros(1, ACsize(2));
TofAC = zeros(1, ACsize(2));
absAC = abs(ACs);
absAC = absAC(2,:);
for i = 1 : ACsize(2)
    if (absAC(i) == 0)
        % Mark the end of the block
        KofAC(i) = -1;
        TofAC(i) = -1;
    else
        KofAC(i) = floor(log2(absAC(i))) + 1;
        TofAC(i) = absAC(i) - pow2(KofAC(i) - 1);
    end
end
% 3. huffman code k
% (1). conbine H and d then calculate the probabilities.
DnK = {};
CountDK = [];
for i = 1 : ACsize(2)
    if ACs(2,i) == 0
    else
        sizeDK = size(DnK);
        sizeDK = sizeDK(2);
        flag = 0;
        if sizeDK == 0
            DnK{sizeDK + 1} = [DofAC(i), KofAC(i)];
            CountDK = [CountDK, 1];
        else
            for k = 1 : sizeDK
                if DnK{k} == [DofAC(i), KofAC(i)];
                    CountDK(k) = CountDK(k) + 1;
                    flag = 1;
                end
            end
            if flag == 0
                DnK{sizeDK + 1} = [DofAC(i), KofAC(i)];
                CountDK = [CountDK, 1];
            end
        end
    end
end
% Add the EOB sign to do Huffman coding
DnK{sizeDK + 1} = 'EOB';
CountDK = [CountDK, numOfBlock];

total = 0;
sizeDK = size(DnK);
sizeDK = sizeDK(2);
for i = 1 : sizeDK
    total = total + CountDK(i);
end
for i = 1 : sizeDK
    CountDK(i) = CountDK(i) / total;
end
% (2). Huffman code
HdictAc = huffmandict(DnK,CountDK);
HofAC1 = cell(1, ACsize(2));
for i = 1 : ACsize(2)
    if ACs(2,i) ~= 0
        tempcell = {[DofAC(i), KofAC(i)]};
        HofAC1{1, i} = huffmanenco(tempcell, HdictAc);
    else
        HofAC1{1, i} = [];
    end
end

% 4. Convert t into binary. 
MofAC = cell(1, ACsize(2));
for i = 1 : ACsize(2)
    if ACs(2,i) == 0   
        continue;
    end
    if KofAC(i) == 1 | KofAC(i) == 0
        temp1 = [];
    else
        temp = de2bi(TofAC(i),KofAC(i)-1);
        [tempRow, tempCol] = size(temp);
        temp1 = [];
        for j = 0 : (tempCol - 1)
            temp1 = [temp1; temp(tempCol - j)];
        end
    end
    MofAC{1, i} = temp1;
end

% ----------     Now AC coefficients can be encoded.     ----------

% ----------     Head of DC coefficients     ----------
% Although in this case, ks that greater than 7 never shows up, I still store them for generalization.
headOfDC = [];
% l is determined by the max length of the huffman code
maxOfHfm = 0;
for i = 1 : 12
    temp = size(Hdict{i,2});
    maxOfHfm = max(maxOfHfm, temp(2));
end

maxOfHfmAC = 0;
for i = 1 : 12
    temp = size(Hdict{i,2});
    maxOfHfmAC = max(maxOfHfmAC, temp(2));
end

maxOfHfm = max(maxOfHfm, maxOfHfmAC);

l=ceil(log2(maxOfHfm));
% convert l to binary in 8 bits.
lbi = de2bi(maxOfHfm, 8);
[lbiRow, lbiCol] = size(lbi);
lbi1 = [];
for j = 0 : (lbiCol - 1)
    lbi1 = [lbi1, lbi(lbiCol - j)];
end
headOfDC = [headOfDC, lbi1];
% Assume there is the agreement: 
% We use 8 bits to represent Symbols.
for i = 1 : 12
    % Symbol
    tempSymbol = de2bi(i , 8);
    [tempSymbolRow, tempSymbolCol] = size(tempSymbol);
    tempSymbol1 = [];
    for j = 0 : (tempSymbolCol - 1)
        tempSymbol1 = [tempSymbol1, tempSymbol(tempSymbolCol - j)];
    end
    % length of code
    lenOflen = size(Hdict{i,2});
    lenOflen = lenOflen(2);
    templen = de2bi(lenOflen, l);
    [templenRow, templenCol] = size(templen);
    templen1 = [];
    for j = 0 : (templenCol - 1)
        templen1 = [templen1, templen(templenCol - j)];
    end
    % code
    tempCode = Hdict{i,2};
    headOfDC = [headOfDC, tempSymbol1, templen1, tempCode];
end


% ----------     Head of AC coefficients     ----------
headOfAC = [];
for i = 1 : 11
    % Symbol
    tempSymbol = de2bi((i + 12) , 8);
    [tempSymbolRow, tempSymbolCol] = size(tempSymbol);
    tempSymbol1 = [];
    for j = 0 : (tempSymbolCol - 1)
        tempSymbol1 = [tempSymbol1, tempSymbol(tempSymbolCol - j)];
    end
    % length of code
    lenOflen = size(HdictAc{i,2});
    lenOflen = lenOflen(2);
    templen = de2bi(lenOflen, l);
    [templenRow, templenCol] = size(templen);
    templen1 = [];
    for j = 0 : (templenCol - 1)
        templen1 = [templen1, templen(templenCol - j)];
    end
    % code
    tempCode = HdictAc{i,2};
    headOfAC = [headOfAC, tempSymbol1, templen1, tempCode];
end

% Bit Stream of DC coefficients
bsDC = [];
for i = 1 : sizeDCdpcm(2)
    bsDC = [bsDC, HofDC1{i}', SofDC(i), MofDC{i}'];
end
% Bit Stream of AC coefficients    
bsAC = [];
for i = 1 : ACsize(2)
    if ACs(2,i) ~= 0
        bsAC = [bsAC, HofAC1{i}', SofAC(i), MofAC{i}'];
    end
end


% ----------     FINAL BIT STREAM     ----------
BS = [headOfDC, headOfAC, bsDC, bsAC];
BSsize = size(BS);
BSsize = BSsize(2);
Hsize = size(HdictAc);
Hsize = Hsize(1);
% The last one in dictionary is EOB accroding to the code.
EOBbit = size(HdictAc{Hsize,2});
% Add bits of EOBs into the total bits manually
CR = (orgnSize + EOBbit(2)* numOfBlock) / BSsize;


% ----------     Partially Decode     ----------
Gde = [];
for i = 1:rows/8
    newRow = [];
    for j = 1:cols/8
        temp = Ground( (i-1)*8+1 : i*8 , (j-1)*8+1 : j*8);
        temp = temp .* (Qparam * Q);
        temp = idct2(temp);
        temp = temp + 128;
        newRow = [newRow, temp];
    end
    Gde = [Gde; newRow];
end
% Test
%imagesc(Gde);
%colormap(gray);

snr(Gde, G)
Qparam
CR