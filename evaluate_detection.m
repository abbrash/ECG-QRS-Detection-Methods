function [sensitivity, positive_predictivity] = evaluate_detection(detected_qrs, annotated_qrs, tolerance)
% Evaluate QRS detection performance
%
% Inputs:
%   detected_qrs: Array of QRS complex locations detected by the algorithm
%   annotated_qrs: Array of true QRS complex locations (ground truth)
%   tolerance: Maximum allowed time difference between detected and true QRS locations (in samples)
%
% Outputs:
%   sensitivity: The proportion of actual QRS complexes that were correctly identified
%   positive_predictivity: The proportion of detected QRS complexes that were true QRS complexes

% Calculate true positives
% A detection is considered correct if it's within the tolerance of an annotation
true_positives = sum(min(abs(detected_qrs(:) - annotated_qrs(:)')) <= tolerance);

% Calculate false negatives (missed detections)
false_negatives = length(annotated_qrs) - true_positives;

% Calculate false positives (spurious detections)
false_positives = length(detected_qrs) - true_positives;

% Calculate sensitivity (also known as recall or true positive rate)
% Sensitivity = TP / (TP + FN)
sensitivity = true_positives / (true_positives + false_negatives);

% Calculate positive predictivity (also known as precision or positive predictive value)
% Positive Predictivity = TP / (TP + FP)
positive_predictivity = true_positives / (true_positives + false_positives);

end