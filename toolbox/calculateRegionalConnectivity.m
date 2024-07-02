function [regionalMatrix] = calculateRegionalConnectivity(FullMatrix, x, EEGchanLocs)

% Read channel locations
chanLocs = cell(length(EEGchanLocs), 1);
for i = 1:length(EEGchanLocs)
    chanLocs{i} = EEGchanLocs(i).labels;
end

% Re-order the matrix based on regions
LeftFrontal_labels = {'Fp1', 'F3', 'F7', 'AF7', 'AF3', 'F1', 'F5'};
LeftFrontal_index = zeros(length(LeftFrontal_labels), 1);
for i = 1:length(LeftFrontal_labels)
    currentLabel = LeftFrontal_labels{i};
    LeftFrontal_index(i) = find(strcmp(currentLabel, chanLocs));
end

RightFrontal_labels = {'F4', 'F8', 'Fp2', 'F6', 'AF8', 'AF4', 'F2'};
RightFrontal_index = zeros(length(RightFrontal_labels), 1);
for i = 1:length(RightFrontal_labels)
    currentLabel = RightFrontal_labels{i};
    RightFrontal_index(i) = find(strcmp(currentLabel, chanLocs));
end

LeftCentral_labels = {'FC5', 'FC1', 'C3', 'CP5', 'CP1', 'FC3', 'C1', 'C5', 'CP3'};
LeftCentral_index = zeros(length(LeftCentral_labels), 1);
for i = 1:length(LeftCentral_labels)
    currentLabel = LeftCentral_labels{i};
    LeftCentral_index(i) = find(strcmp(currentLabel, chanLocs));
end

RightCentral_labels = {'CP6', 'CP2', 'C4', 'FC6', 'FC2', 'CP4', 'C6', 'C2', 'FC4'};
RightCentral_index = zeros(length(RightCentral_labels), 1);
for i = 1:length(RightCentral_labels)
    currentLabel = RightCentral_labels{i};
    RightCentral_index(i) = find(strcmp(currentLabel, chanLocs));
end

LeftTemporal_labels = {'FT9', 'T7', 'TP9', 'P7', 'FT7', 'TP7'};
LeftTemporal_index = zeros(length(LeftTemporal_labels), 1);
for i = 1:length(LeftTemporal_labels)
    currentLabel = LeftTemporal_labels{i};
    LeftTemporal_index(i) = find(strcmp(currentLabel, chanLocs));
end

RightTemporal_labels = {'P8', 'TP10', 'T8', 'FT10', 'TP8', 'FT8'};
RightTemporal_index = zeros(length(RightTemporal_labels), 1);
for i = 1:length(RightTemporal_labels)
    currentLabel = RightTemporal_labels{i};
    RightTemporal_index(i) = find(strcmp(currentLabel, chanLocs));
end

LeftParieto_labels = {'P3', 'O1', 'P1', 'P5', 'PO7', 'PO3'};
LeftParieto_index = zeros(length(LeftParieto_labels), 1);
for i = 1:length(LeftParieto_labels)
    currentLabel = LeftParieto_labels{i};
    LeftParieto_index(i) = find(strcmp(currentLabel, chanLocs));
end

RightParieto_labels = {'O2', 'P4', 'PO4', 'PO8', 'P6', 'P2'};
RightParieto_index = zeros(length(RightParieto_labels), 1);
for i = 1:length(RightParieto_labels)
    currentLabel = RightParieto_labels{i};
    RightParieto_index(i) = find(strcmp(currentLabel, chanLocs));
end

if ~isempty(find(strcmp('Fz', chanLocs), 1))
    Vertex_labels = {'Pz', 'Oz', 'Fz', 'AFz', 'POz', 'CPz', 'FCz'}; 
else
    Vertex_labels = {'Pz', 'Oz', 'Cz', 'AFz', 'POz', 'CPz', 'FCz'}; 
end
Vertex_index = zeros(length(Vertex_labels), 1);
for i = 1:length(Vertex_labels)
    currentLabel = Vertex_labels{i};
    Vertex_index(i) = find(strcmp(currentLabel, chanLocs));
end

%Take a 63*63 Matrix and change to 9*9 Matrix
LFLF = mean(FullMatrix(LeftFrontal_index, LeftFrontal_index),[1,2]);
LFRF = mean(FullMatrix(LeftFrontal_index, RightFrontal_index),[1,2]);
LFLC = mean(FullMatrix(LeftFrontal_index, LeftCentral_index),[1,2]);
LFRC = mean(FullMatrix(LeftFrontal_index, RightCentral_index),[1,2]);
LFLT = mean(FullMatrix(LeftFrontal_index, LeftTemporal_index),[1,2]);
LFRT = mean(FullMatrix(LeftFrontal_index, RightTemporal_index),[1,2]);
LFLP = mean(FullMatrix(LeftFrontal_index, LeftParieto_index),[1,2]);
LFRP = mean(FullMatrix(LeftFrontal_index, RightParieto_index),[1,2]);
LFV = mean(FullMatrix(LeftFrontal_index, Vertex_index),[1,2]);

RFLF = mean(FullMatrix(RightFrontal_index, LeftFrontal_index),[1,2]);
RFRF = mean(FullMatrix(RightFrontal_index, RightFrontal_index),[1,2]);
RFLC = mean(FullMatrix(RightFrontal_index, LeftCentral_index),[1,2]);
RFRC = mean(FullMatrix(RightFrontal_index, RightCentral_index),[1,2]);
RFLT = mean(FullMatrix(RightFrontal_index, LeftTemporal_index),[1,2]);
RFRT = mean(FullMatrix(RightFrontal_index, RightTemporal_index),[1,2]);
RFLP = mean(FullMatrix(RightFrontal_index, LeftParieto_index),[1,2]);
RFRP = mean(FullMatrix(RightFrontal_index, RightParieto_index),[1,2]);
RFV = mean(FullMatrix(RightFrontal_index, Vertex_index),[1,2]);

LCLF = mean(FullMatrix(LeftCentral_index, LeftFrontal_index),[1,2]);
LCRF = mean(FullMatrix(LeftCentral_index, RightFrontal_index),[1,2]);
LCLC = mean(FullMatrix(LeftCentral_index, LeftCentral_index),[1,2]);
LCRC = mean(FullMatrix(LeftCentral_index, RightCentral_index),[1,2]);
LCLT = mean(FullMatrix(LeftCentral_index, LeftTemporal_index),[1,2]);
LCRT = mean(FullMatrix(LeftCentral_index, RightTemporal_index),[1,2]);
LCLP = mean(FullMatrix(LeftCentral_index, LeftParieto_index),[1,2]);
LCRP = mean(FullMatrix(LeftCentral_index, RightParieto_index),[1,2]);
LCV = mean(FullMatrix(LeftCentral_index, Vertex_index),[1,2]);

RCLF = mean(FullMatrix(RightCentral_index, LeftFrontal_index),[1,2]);
RCRF = mean(FullMatrix(RightCentral_index, RightFrontal_index),[1,2]);
RCLC = mean(FullMatrix(RightCentral_index, LeftCentral_index),[1,2]);
RCRC = mean(FullMatrix(RightCentral_index, RightCentral_index),[1,2]);
RCLT = mean(FullMatrix(RightCentral_index, LeftTemporal_index),[1,2]);
RCRT = mean(FullMatrix(RightCentral_index, RightTemporal_index),[1,2]);
RCLP = mean(FullMatrix(RightCentral_index, LeftParieto_index),[1,2]);
RCRP = mean(FullMatrix(RightCentral_index, RightParieto_index),[1,2]);
RCV = mean(FullMatrix(RightCentral_index, Vertex_index),[1,2]);

LTLF = mean(FullMatrix(LeftTemporal_index, LeftFrontal_index),[1,2]);
LTRF = mean(FullMatrix(LeftTemporal_index, RightFrontal_index),[1,2]);
LTLC = mean(FullMatrix(LeftTemporal_index, LeftCentral_index),[1,2]);
LTRC = mean(FullMatrix(LeftTemporal_index, RightCentral_index),[1,2]);
LTLT = mean(FullMatrix(LeftTemporal_index, LeftTemporal_index),[1,2]);
LTRT = mean(FullMatrix(LeftTemporal_index, RightTemporal_index),[1,2]);
LTLP = mean(FullMatrix(LeftTemporal_index, LeftParieto_index),[1,2]);
LTRP = mean(FullMatrix(LeftTemporal_index, RightParieto_index),[1,2]);
LTV = mean(FullMatrix(LeftTemporal_index, Vertex_index),[1,2]);

RTLF = mean(FullMatrix(RightTemporal_index, LeftFrontal_index),[1,2]);
RTRF = mean(FullMatrix(RightTemporal_index, RightFrontal_index),[1,2]);
RTLC = mean(FullMatrix(RightTemporal_index, LeftCentral_index),[1,2]);
RTRC = mean(FullMatrix(RightTemporal_index, RightCentral_index),[1,2]);
RTLT = mean(FullMatrix(RightTemporal_index, LeftTemporal_index),[1,2]);
RTRT = mean(FullMatrix(RightTemporal_index, RightTemporal_index),[1,2]);
RTLP = mean(FullMatrix(RightTemporal_index, LeftParieto_index),[1,2]);
RTRP = mean(FullMatrix(RightTemporal_index, RightParieto_index),[1,2]);
RTV = mean(FullMatrix(RightTemporal_index, Vertex_index),[1,2]);

LPLF = mean(FullMatrix(LeftParieto_index, LeftFrontal_index),[1,2]);
LPRF = mean(FullMatrix(LeftParieto_index, RightFrontal_index),[1,2]);
LPLC = mean(FullMatrix(LeftParieto_index, LeftCentral_index),[1,2]);
LPRC = mean(FullMatrix(LeftParieto_index, RightCentral_index),[1,2]);
LPLT = mean(FullMatrix(LeftParieto_index, LeftTemporal_index),[1,2]);
LPRT = mean(FullMatrix(LeftParieto_index, RightTemporal_index),[1,2]);
LPLP = mean(FullMatrix(LeftParieto_index, LeftParieto_index),[1,2]);
LPRP = mean(FullMatrix(LeftParieto_index, RightParieto_index),[1,2]);
LPV = mean(FullMatrix(LeftParieto_index, Vertex_index),[1,2]);

RPLF = mean(FullMatrix(RightParieto_index, LeftFrontal_index),[1,2]);
RPRF = mean(FullMatrix(RightParieto_index, RightFrontal_index),[1,2]);
RPLC = mean(FullMatrix(RightParieto_index, LeftCentral_index),[1,2]);
RPRC = mean(FullMatrix(RightParieto_index, RightCentral_index),[1,2]);
RPLT = mean(FullMatrix(RightParieto_index, LeftTemporal_index),[1,2]);
RPRT = mean(FullMatrix(RightParieto_index, RightTemporal_index),[1,2]);
RPLP = mean(FullMatrix(RightParieto_index, LeftParieto_index),[1,2]);
RPRP = mean(FullMatrix(RightParieto_index, RightParieto_index),[1,2]);
RPV = mean(FullMatrix(RightParieto_index, Vertex_index),[1,2]);

VLF = mean(FullMatrix(Vertex_index, LeftFrontal_index),[1,2]);
VRF = mean(FullMatrix(Vertex_index, RightFrontal_index),[1,2]);
VLC = mean(FullMatrix(Vertex_index, LeftCentral_index),[1,2]);
VRC = mean(FullMatrix(Vertex_index, RightCentral_index),[1,2]);
VLT = mean(FullMatrix(Vertex_index, LeftTemporal_index),[1,2]);
VRT = mean(FullMatrix(Vertex_index, RightTemporal_index),[1,2]);
VLP = mean(FullMatrix(Vertex_index, LeftParieto_index),[1,2]);
VRP = mean(FullMatrix(Vertex_index, RightParieto_index),[1,2]);
VV = mean(FullMatrix(Vertex_index, Vertex_index),[1,2]);

if x==1 %1=left, 2=right

    %     regionalMatrix = [LFLF LFRF LFLC LFRC LFLT LFRT LFLP LFRP LFV;
    %         RFLF RFRF RFLC RFRC RFLT RFRT RFLP RFRP RFV;
    %         LCLF LCRF LCLC LCRC LCLT LCRT LCLP LCRP LCV;
    %         RCLF RCRF RCLC RCRC RCLT RCRT RCLP RCRP RCV;
    %         LTLF LTRF LTLC LTRC LTLT LTRT LTLP LTRP LTV;
    %         RTLF RTRF RTLC RTRC RTLT RTRT RTLP RTRP RTV;
    %         LPLF LPRF LPLC LPRC LPLT LPRT LPLP LPRP LPV;
    %         RPLF RPRF RPLC RPRC RPLT RPRT RPLP RPRP RPV;
    %         VLF VRF VLC VRC VLT VRT VLP VRP VV];

    % Drop Vertex
    regionalMatrix = [LFLF LFRF LFLC LFRC LFLT LFRT LFLP LFRP;
        RFLF RFRF RFLC RFRC RFLT RFRT RFLP RFRP;
        LCLF LCRF LCLC LCRC LCLT LCRT LCLP LCRP;
        RCLF RCRF RCLC RCRC RCLT RCRT RCLP RCRP;
        LTLF LTRF LTLC LTRC LTLT LTRT LTLP LTRP;
        RTLF RTRF RTLC RTRC RTLT RTRT RTLP RTRP;
        LPLF LPRF LPLC LPRC LPLT LPRT LPLP LPRP;
        RPLF RPRF RPLC RPRC RPLT RPRT RPLP RPRP];

elseif x==2

    %     regionalMatrix = [RFRF RFLF RFRC RFLC RFRT RFLT RFRP RFLP RFV;
    %         LFRF LFLF LFRC LFLC LFRT LFLT LFRP LFLP LFV;
    %         RCRF RCLF RCRC RCLC RCRT RCLT RCRP RCLP RCV;
    %         LCRF LCLF LCRC LCLC LCRT LCLT LCRP LCLP LCV;
    %         RTRF RTLF RTRC RTLC RTRT RTLT RTRP RTLP RTV;
    %         LTRF LTLF LTRC LTLC LTRT LTLT LTRP LTLP LTV;
    %         RPRF RPLF RPRC RPLC RPRT RPLT RPRP RPLP RPV;
    %         LPRF LPLF LPRC LPLC LPRT LPLT LPRP LPLP LPV;
    %         VRF VLF VRC VLC VRT VLT VRP VLP VV];

    % Drop Vertex
    regionalMatrix = [RFRF RFLF RFRC RFLC RFRT RFLT RFRP RFLP;
        LFRF LFLF LFRC LFLC LFRT LFLT LFRP LFLP;
        RCRF RCLF RCRC RCLC RCRT RCLT RCRP RCLP;
        LCRF LCLF LCRC LCLC LCRT LCLT LCRP LCLP;
        RTRF RTLF RTRC RTLC RTRT RTLT RTRP RTLP;
        LTRF LTLF LTRC LTLC LTRT LTLT LTRP LTLP;
        RPRF RPLF RPRC RPLC RPRT RPLT RPRP RPLP;
        LPRF LPLF LPRC LPLC LPRT LPLT LPRP LPLP];

end

