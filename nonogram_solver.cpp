#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
//using namespace std;
/*--------------------------------------------------------------------------------//Définitions//--------------------------------------------------------------------------------*/
#define INDETERMINE 0
#define BLACK100 1
#define WHITE100 2
#define BLACKMAYBE 3
#define MISMATCHPENALTY -1
#define GAPPENALTY -1
#define ERRORACCEPTANCEREDUCTOR 1.1
float errorAcceptance = 0.1; // 1 = accépté 100% d erreur

/*--------------------------------------------------------------------------------//CODE//--------------------------------------------------------------------------------*/
int RowSize = 0;
int ColSize = 0;
std::string filename = "pti\\10x10.pti"; // Chemin du fichier

void printMatrix(const std::vector<std::vector<int>>& matrix) {
    for (const auto& row : matrix) {
        for (int val : row) {
            std::cout << val << "\t";
        }
        std::cout << std::endl;
    }
}

void needlemanWunsch(std::vector<short>& sequence1, std::vector<short>& sequence2) {

    int m = sequence1.size();
    int n = sequence2.size();

    // Initialize the dynamic programming matrix
    std::vector<std::vector<int>> dpMatrix(m + 1, std::vector<int>(n + 1, 0));

    // Initialize the first row and column with gap penalties
    for (int i = 1; i <= m; ++i) {
        dpMatrix[i][0] = dpMatrix[i - 1][0] + GAPPENALTY;
    }

    for (int j = 1; j <= n; ++j) {
        dpMatrix[0][j] = dpMatrix[0][j - 1] + GAPPENALTY;
    }

    // Fill the dynamic programming matrix
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int match = dpMatrix[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? sequence1[i - 1] * sequence1[i - 1] : MISMATCHPENALTY);
            int gap1 = dpMatrix[i - 1][j] + GAPPENALTY;
            int gap2 = dpMatrix[i][j - 1] + GAPPENALTY;

            dpMatrix[i][j] = std::max({ match, gap1, gap2 });
        }
    }
    /*
    // Print the dynamic programming matrix
    std::cout << "Dynamic Programming Matrix:" << std::endl;
    printMatrix(dpMatrix);
    std::cout << std::endl;
    */


    // Backtrack to find the optimal alignment
    int i = m, j = n;
    std::vector<short> alignedSequence1, alignedSequence2;
    //vector

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && dpMatrix[i][j] == dpMatrix[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? sequence1[i - 1] * sequence1[i - 1] : MISMATCHPENALTY)) {
            alignedSequence1.insert(alignedSequence1.begin(), sequence1[i - 1]);
            alignedSequence2.insert(alignedSequence2.begin(), sequence2[j - 1]);
            --i;
            --j;
        }
        else if (i > 0 && dpMatrix[i][j] == dpMatrix[i - 1][j] + GAPPENALTY) {
            alignedSequence1.insert(alignedSequence1.begin(), sequence1[i - 1]);
            alignedSequence2.insert(alignedSequence2.begin(), 0);
            --i;
        }
        else {
            alignedSequence1.insert(alignedSequence1.begin(), 0);
            alignedSequence2.insert(alignedSequence2.begin(), sequence2[j - 1]);
            --j;
        }
    }
    sequence1 = alignedSequence1;
    sequence2 = alignedSequence2;
    /*
    // Print the optimal alignment
    std::cout << std::endl << "Optimal Alignment:" << std::endl;
    for (int num : alignedSequence1) std::cout << num << " ";
    std::cout << std::endl;
    for (int num : alignedSequence2) std::cout << num << " ";
    std::cout << std::endl;
    */
}

bool movingRight(std::vector<short>& tempRowEnd, std::vector<short>& tempRowStart, std::vector<short>& MovedPlace, short seedGrp)
{
    if (seedGrp < tempRowStart.size() - 1) // si pas dernier GRP
    {
        if (tempRowStart[seedGrp + 1] - tempRowEnd[seedGrp] > 2) // SI espace de bougé
        {
            for (short i = tempRowStart[seedGrp]; i < tempRowEnd[seedGrp] + 2; i++)// Endroits changé
            {
                MovedPlace.push_back(i);
                //std::cout << std::endl << "MovedPlace :" << i << std::endl;
            }
            //group can move
            tempRowStart[seedGrp] += 1;
            tempRowEnd[seedGrp] += 1;
            return true;
        }
        else { // Suivant trop proche
            // MOVING NEXT FIRST ?
            if (movingRight(tempRowEnd, tempRowStart, MovedPlace, seedGrp + 1)) {
                for (short i = tempRowStart[seedGrp]; i < tempRowEnd[seedGrp] + 2; i++)// Endroits changé
                {
                    MovedPlace.push_back(i);
                    //std::cout << std::endl << "MovedPlace :" << i << std::endl;

                }
                //group can move
                tempRowStart[seedGrp] += 1;
                tempRowEnd[seedGrp] += 1;
                return true;
            }
            else {
                return false;
            }
        }
    }
    else {// si Dernier GRP
        if (tempRowEnd[seedGrp] < (RowSize - 1))
        {
            for (short i = tempRowStart[seedGrp]; i < tempRowEnd[seedGrp] + 2; i++)// Endroits changé
            {
                MovedPlace.push_back(i);
                //std::cout << std::endl << "MovedPlace :"<< i << std::endl;

            }
            //group can move
            tempRowStart[seedGrp] += 1;
            tempRowEnd[seedGrp] += 1;
            return true;
        }
        else {// JUMP ? OR NO MOVE ? 
            //NO MOVE
            return false;
        }
    }
}

bool movingLeft(std::vector<short>& tempRowEnd, std::vector<short>& tempRowStart, std::vector<short>& MovedPlace, short seedGrp)
{
    if (seedGrp > 0) // si pas premier GRP
    {
        if (tempRowStart[seedGrp] - tempRowEnd[seedGrp - 1] > 2) // SI espace de bougé
        {
            for (short i = tempRowStart[seedGrp] - 1; i < tempRowEnd[seedGrp] + 1; i++)// Endroits changé
            {
                MovedPlace.push_back(i);
                //std::cout << std::endl << "MovedPlace :" << i << std::endl;
            }
            //group can move
            tempRowStart[seedGrp] -= 1;
            tempRowEnd[seedGrp] -= 1;
            return true;
        }
        else { // Suivant trop proche
            // MOVING NEXT FIRST ?
            if (movingLeft(tempRowEnd, tempRowStart, MovedPlace, seedGrp - 1)) {
                for (short i = tempRowStart[seedGrp] - 1; i < tempRowEnd[seedGrp] + 1; i++)// Endroits changé
                {
                    MovedPlace.push_back(i);
                    //std::cout << std::endl << "MovedPlace :" << i << std::endl;
                }
                //group can move
                tempRowStart[seedGrp] -= 1;
                tempRowEnd[seedGrp] -= 1;
                return true;
            }
            else {
                return false;
            }
        }
    }
    else {// si Premier GRP
        if (tempRowStart[seedGrp] > 0)
        {
            for (short i = tempRowStart[seedGrp] - 1; i < tempRowEnd[seedGrp] + 1; i++)// Endroits changé
            {
                MovedPlace.push_back(i);
                //std::cout << std::endl << "MovedPlace :" << i << std::endl;
            }
            //group can move
            tempRowStart[seedGrp] -= 1;
            tempRowEnd[seedGrp] -= 1;
            return true;
            return true;
        }
        else {// JUMP ? OR NO MOVE ? 
            //NO MOVE
            return false;
        }
    }
}

bool fitModif(short* tableauOrig, short* tableauProp)
{
    for (short j = 0; j < ColSize; j++)
    {
        if (tableauProp[j] == BLACKMAYBE) {
            if (tableauOrig[j] == WHITE100) {
                return false;
            }
        }
        else {
            if (tableauOrig[j] == BLACK100) {
                return false;
            }
        }
    }
    return true;
}

void MoveGroupsFirst(short* line, std::vector<short>& contraintesRows) {
    std::vector<short> gStart;
    // Detection des groupements
    if (line[0] == BLACKMAYBE) {
        // start new grp de 1 a la fin 
        gStart.push_back(0);
    }
    for (short i = 1; i < RowSize; ++i) {
        if (line[i] == BLACKMAYBE && line[i - 1] == INDETERMINE) {
            // start new grp
            gStart.push_back(i);

        }
    }
    std::vector<short> gEnd;
    for (short i = 0; i < gStart.size(); ++i) {
        gEnd.push_back(gStart[i] + contraintesRows[i] - 1);
    }
    // Vérifie si le dernier groupement peut être déplacé
    short cursor = gStart.size() - 1;
    if (gEnd[cursor] < RowSize - 1)
    {
        //last group can move
        gStart[cursor] += 1;
        gEnd[cursor] += 1;
    }
    else {
        cursor -= 1;
        bool notDecal = true;
        while (notDecal) {
            if (gStart[cursor + 1] - gEnd[cursor] > 2) {
                //movable
                gStart[cursor] += 1;
                gEnd[cursor] += 1;
                for (short i = cursor + 1; i < gStart.size(); ++i) {
                    gStart[i] = gEnd[i - 1] + 2;
                    gEnd[i] = gStart[i] + contraintesRows[i] - 1;
                }
                notDecal = false;
            }
            else {
                cursor -= 1;
                //retry
            }
        }
    }

    for (short j = 0; j < RowSize; ++j) {
        line[j] = INDETERMINE;
    }
    for (short i = 0; i < gStart.size(); ++i) {
        for (short j = gStart[i]; j < gEnd[i] + 1; ++j) {
            line[j] = BLACKMAYBE;
        }
    }
}

void createPTOFile(const char* InFilePath, short** matrix, int score, int rows, int cols) {
    std::string filePath = "ptoOutput\\" + std::string(InFilePath) + ".pto";
    std::ofstream outFile(filePath);
    if (!outFile) {
        std::cerr << "Erreur lors de la création du fichier PTO." << std::endl;
        return;
    }

    outFile << "Name: Test Solution " << std::string(InFilePath) << "\n";
    outFile << "Solver: TestSolver\n";
    outFile << "Rows: " << rows << "\n";
    outFile << "Cols: " << cols << "\n";
    outFile << "Score: " << score << "\n\n";

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Appliquer la transformation
            short valueToWrite;
            if (matrix[i][j] == 3) {
                valueToWrite = 1;
            }
            else if (matrix[i][j] == 2) {
                valueToWrite = 0;
            }
            else {
                valueToWrite = matrix[i][j];
            }

            outFile << valueToWrite;
        }
        outFile << "\n";
    }
}

int main()
{

    std::vector<std::vector<short>> contraintesCols;
    std::vector<std::vector<short>> contraintesRows;
    std::string cuttedName = filename.substr(filename.find_last_of("\\/") + 1); // Nom avec extension
    cuttedName = cuttedName.substr(0, cuttedName.find_last_of(".")); // Retirer l'extension
    const char* cuttedNamechar = cuttedName.c_str();
    std::fstream myfile(filename);  // Remplacez par le chemin correct
    std::string myline;
    short nbrenoir = 0;
    bool firstEmptyLine = false;  // Flag pour la première ligne vide
    bool secondEmptyLine = false;  // Flag pour la deuxième ligne vide

    while (getline(myfile, myline)) {
        if (myline.rfind("Rows:", 0) == 0) {
            sscanf_s(myline.c_str(), "Rows: %d", &RowSize);
        }
        if (myline.rfind("Cols:", 0) == 0) {
            sscanf_s(myline.c_str(), "Cols: %d", &ColSize);
        }

        if (myline.empty()) {
            if (!firstEmptyLine) {
                firstEmptyLine = true;  // On a rencontré la première ligne vide
                continue;
            }
            else if (!secondEmptyLine) {
                secondEmptyLine = true;  // On a rencontré la deuxième ligne vide
                continue;
            }
            else {
                break;  // On a déjà passé les deux lignes vides, on sort
            }
        }

        if (firstEmptyLine && !secondEmptyLine) {
            std::stringstream ss(myline);
            std::vector<short> constraints;
            short value;
            while (ss >> value) {
                nbrenoir += value;
                constraints.push_back(value);
            }
            contraintesCols.push_back(constraints);  // Ajouter aux colonnes
        }

        if (secondEmptyLine) {
            std::stringstream ss(myline);
            std::vector<short> constraints;
            short value;
            while (ss >> value) {
                constraints.push_back(value);
            }
            contraintesRows.push_back(constraints);  // Ajouter aux lignes
        }
    }

    myfile.close();

    //DynamiqueVersion
    short** tableau2D = new short* [RowSize];
    for (short i = 0; i < RowSize; i++) {
        tableau2D[i] = new short[ColSize];
        for (short j = 0; j < ColSize; j++) {
            tableau2D[i][j] = INDETERMINE;
        }
    }

    /*--------------------------------------------------------------------------------//Pre-Process//--------------------------------------------------------------------------------*/

    //ROW PREPROCESS
    short HalfRows = RowSize / 2;// Hard Rounded
    for (short i = 0; i < RowSize; i++)
    {
        //count whitespaces
        short tmp = contraintesRows[i].size() - 1;

        for (short j = 0; j < contraintesRows[i].size(); j++)
        {
            //count black groupe
            tmp += contraintesRows[i][j];
        }
        // ATTENTION ELSE IF FULL -> faster and also white  ????????????????????????
        // Si sommme plus grande que la moitié ->  MAYBE plassable
        if (tmp > HalfRows)
        {
            // si contraintesRows > (Cols - tmp) -> plassable
            short VariableSpace = 0;
            for (short j = 0; j < contraintesRows[i].size(); j++)
            {
                for (short k = 0; k < contraintesRows[i][j]; k++)
                {
                    if (k >= (RowSize - tmp))
                    {
                        tableau2D[i][k + VariableSpace] = BLACK100;
                    }
                }
                VariableSpace += contraintesRows[i][j] + 1; // +1 for Whitespace 
            }
        }
    }

    //Cols PREPROCESS
    short HalfCols = ColSize / 2;// Hard Rounded
    for (short i = 0; i < ColSize; i++)
    {
        //count whitespaces
        short tmp = contraintesCols[i].size() - 1;

        for (short j = 0; j < contraintesCols[i].size(); j++)
        {
            //count black groupe
            tmp += contraintesCols[i][j];
        }
        // ATTENTION ELSE IF FULL -> faster and also white  ????????????????????????
        // Si sommme plus grande que la moitié ->  MAYBE plassable
        if (tmp > HalfCols)
        {
            // si contraintesRows > (Cols - tmp) -> plassable
            short VariableSpace = 0;
            for (short j = 0; j < contraintesCols[i].size(); j++)
            {
                for (short k = 0; k < contraintesCols[i][j]; k++)
                {
                    if (k >= (ColSize - tmp))
                    {
                        tableau2D[k + VariableSpace][i] = BLACK100;
                    }
                }
                VariableSpace += contraintesCols[i][j] + 1; // +1 for Whitespace 
            }
        }
    }



    // AFFICHAGE SOL
    std::cout << std::endl;
    std::cout << "AFTER PREPROCESS";
    for (short i = 0; i < RowSize; i++)
    {
        std::cout << std::endl;
        for (short j = 0; j < ColSize; j++)
        {
            std::cout << tableau2D[i][j];
        }
    }


    /*--------------------------------------------------------------------------------//Placement des groupements//--------------------------------------------------------------------------------*/
    // rajouter un if full pour pas tourné pour rien ????????????????????????
    std::vector<std::vector<short>> gStart;
    std::vector<std::vector<short>> gEnd;
    gStart.resize(ColSize);
    gEnd.resize(ColSize);
    // pour chaque ligne
    for (short i = 0; i < RowSize; i++)
    {
        // créé un tableau 1D avec les groupements
        short* tableau1DG = new short[RowSize];
        for (short j = 0; j < RowSize; j++) {
            tableau1DG[j] = INDETERMINE;
        }
        short VariableSpace = 0;
        for (short j = 0; j < contraintesRows[i].size(); j++)
        {
            // START GROUPEMENT
            gStart[i].push_back(VariableSpace);
            for (short k = 0; k < contraintesRows[i][j]; k++)
            {
                tableau1DG[k + VariableSpace] = BLACKMAYBE;
            }

            // END GROUPEMENT
            gEnd[i].push_back(contraintesRows[i][j] + VariableSpace - 1);
            VariableSpace += contraintesRows[i][j] + 1; // +1 for Whitespace 
        }
        if (fitModif(tableau2D[i], tableau1DG)) {
            // Appropriation
            for (short j = 0; j < ColSize; j++)
            {
                if (tableau2D[i][j] != BLACK100 && tableau2D[i][j] != WHITE100) {
                    tableau2D[i][j] = tableau1DG[j];
                }
            }
        }
        else {
            //refaire une proposition
            bool notOK = true;

            while (notOK) {
                // decalage des groupements
                // Vérifie si le dernier groupement peut être déplacé
                short cursor = gStart[i].size() - 1;
                if (gEnd[i][cursor] < RowSize - 1)
                {
                    //last group can move
                    gStart[i][cursor] += 1;
                    gEnd[i][cursor] += 1;
                }
                else {
                    cursor -= 1;
                    bool notDecal = true;
                    while (notDecal) {
                        if (gStart[i][cursor + 1] - gEnd[i][cursor] > 2) {
                            //movable
                            gStart[i][cursor] += 1;
                            gEnd[i][cursor] += 1;
                            for (short j = cursor + 1; j < gStart[i].size(); ++j) {
                                gStart[i][j] = gEnd[i][j - 1] + 2;
                                gEnd[i][j] = gStart[i][j] + contraintesRows[i][j] - 1;
                            }
                            notDecal = false;
                        }
                        else {
                            cursor -= 1;
                            //retry
                        }
                    }
                }
                for (short j = 0; j < RowSize; ++j) {
                    tableau1DG[j] = INDETERMINE;
                }
                for (short j = 0; j < gStart[i].size(); ++j) {
                    for (short k = gStart[i][j]; k < gEnd[i][j] + 1; ++k) {
                        tableau1DG[k] = BLACKMAYBE;
                    }
                }
                if (fitModif(tableau2D[i], tableau1DG)) {
                    // Appropriation
                    for (short j = 0; j < ColSize; j++)
                    {
                        if (tableau2D[i][j] != BLACK100 && tableau2D[i][j] != WHITE100) {
                            tableau2D[i][j] = tableau1DG[j];
                        }
                    }
                    notOK = false;
                }
            }
        }
        delete[] tableau1DG;
    }
    // AFFICHAGE Placement
    std::cout << std::endl << std::endl;
    std::cout << "Placement initial";
    for (short i = 0; i < RowSize; i++)
    {
        std::cout << std::endl;
        for (short j = 0; j < ColSize; j++)
        {
            std::cout << tableau2D[i][j];
        }
    }

    /*--------------------------------------------------------------------------------//Recuit Simulé//--------------------------------------------------------------------------------*/
    //// calculate starting Scoring ? 
    // Detection des groupements
    std::vector<short> ScoresCols;
    ScoresCols.resize(RowSize, 0);
    std::vector<short> gdetect;
    bool onTrack = false;
    short TrackLen = 0;
    int ScoreTotal = 0;
    for (short i = 0; i < RowSize; i++) {
        gdetect.clear();
        for (short j = 0; j < ColSize; j++) {
            if (tableau2D[j][i] == WHITE100 || tableau2D[j][i] == INDETERMINE) {//si j arrive sur un blanc
                if (onTrack) {//si je suis(SUIVRE) un noir
                    // FIN DE GROUPEMENT
                    onTrack = false;
                    gdetect.push_back(TrackLen);
                    TrackLen = 0;
                }
            }
            else {
                //pas blanc donc je suis noir
                onTrack = true;
                TrackLen++;
            }
        }
        if (TrackLen != 0) {
            onTrack = false;
            gdetect.push_back(TrackLen);
            TrackLen = 0;
        }
        /*
        // comparer gdetect a connstraints
        std::cout << std::endl << std::endl;
        for (short j = 0; j < gdetect.size(); j++) {
            std::cout << gdetect[j] << " ";
        }
        std::cout << std::endl;
        for (short j = 0; j < contraintesCols[i].size(); j++) {
            std::cout << contraintesCols[i][j] << " ";
        }
        */

        // Version Needle 
        signed short deltaNW = contraintesCols[i].size() - gdetect.size();
        //////////// DEBUG ////////////
        std::cout << std::endl << "delta : " << deltaNW;
        //////////// DEBUG END ////////////

        std::vector<short> seq1 = contraintesCols[i];
        std::vector<short> seq2 = gdetect;
        needlemanWunsch(seq1, seq2);
        for (short j = 0; j < seq1.size(); j++) {
            ScoresCols[i] += pow((seq1[j] - seq2[j]), 2);
        }
        ScoresCols[i] += pow(deltaNW, 2);
        ScoreTotal += ScoresCols[i];
        //////////// DEBUG ////////////
        std::cout << std::endl;
        std::cout << "score : " << ScoresCols[i];
        //////////// DEBUG END ////////////
    }
    //////////// DEBUG ////////////
    std::cout << std::endl;
    std::cout << "scoreTotal : " << ScoreTotal;
    //////////// DEBUG END ////////////

    //// Moving Neighboor Loop

    // Initialise le générateur avec une graine aléatoire
    std::random_device rd; // Fournit une graine (source d'entropie)
    std::mt19937 gen(rd()); // Générateur Mersenne Twister

    std::vector<short> tempRowStart;
    std::vector<short> tempRowEnd;
    std::vector<short> MovedPlace;

    //DynamiqueVersion
    short** tableau2DReplica = new short* [RowSize];
    for (short i = 0; i < RowSize; i++) {
        tableau2DReplica[i] = new short[ColSize];
        for (short j = 0; j < ColSize; j++) {
            tableau2DReplica[i][j] = tableau2D[i][j];
        }
    }


    bool RecuitRunStep = true;
    int iteration = 0;
    int TrueIteration = 0;
    while (RecuitRunStep) {
        std::cout << std::endl << std::endl << "---// Iteration " << iteration << " //--- ";
        MovedPlace.clear();
        //Seed Move Generation
        std::uniform_int_distribution<> distribRow(0, ColSize - 1);
        short seedRow = distribRow(gen);
        //seedRow = 9;
        std::cout << std::endl << std::endl << "SeedRow " << seedRow << " ";
        std::uniform_int_distribution<> distribGrp(0, contraintesRows[seedRow].size() - 1);
        short seedGrp = distribGrp(gen);
        //seedGrp = 1;
        std::cout << std::endl << "SeedGRP " << seedGrp << " ";
        std::uniform_int_distribution<> distribDir(0, 9);
        short seedDir = distribDir(gen);
        //seedDir = 1;
        std::cout << std::endl << "SeedDir " << seedDir << " ";


        // decalge OU PAS 
        tempRowStart = gStart[seedRow];
        tempRowEnd = gEnd[seedRow];
        bool moved = false;
        if (seedDir > 3) { // DROITE
            moved = movingRight(tempRowEnd, tempRowStart, MovedPlace, seedGrp);
        }
        else { // GAUCHE
            moved = movingLeft(tempRowEnd, tempRowStart, MovedPlace, seedGrp);
        }
        if (moved) { // random seed applicable
            std::sort(MovedPlace.begin(), MovedPlace.end());
            // Applicate Change to Replica 
            for (short j = 0; j < RowSize; ++j) {
                tableau2DReplica[seedRow][j] = INDETERMINE;
            }
            for (short j = 0; j < tempRowStart.size(); ++j) {
                for (short k = tempRowStart[j]; k < tempRowEnd[j] + 1; ++k) {

                    tableau2DReplica[seedRow][k] = BLACKMAYBE;
                }
            }
            /*
            // AFFICHAGE SOL
            std::cout << std::endl << std::endl;
            std::cout << "TempState";
            for (short i = 0; i < RowSize; i++)
            {
                std::cout << std::endl;
                for (short j = 0; j < ColSize; j++)
                {
                    std::cout << tableau2DReplica[i][j];
                }
            }
            */


            // Changement applicable ? 
            bool modifjouable = true;
            for (short k = 0; k < MovedPlace.size(); ++k) {

                if (!fitModif(tableau2D[seedRow], tableau2DReplica[seedRow])) {
                    modifjouable = false;
                    std::cout << std::endl << "ERROR MATCHING Bad seed";
                    goto endFitModif;
                }
            }
        endFitModif:
            if (modifjouable) {
                int ActualScore = 0;
                std::vector<short> ScoresModified;
                ScoresModified.resize(MovedPlace.size(), 0);
                int TOTScoresModified = 0;
                //////////// DEBUG ////////////
                std::cout << std::endl << "MATCHING Movable Good seed";
                //////////// DEBUG END ////////////
                // Calculate changes Scores
                //// calculate starting Scoring ? 
                     // Detection des groupements

                std::vector<short> gdetect2;
                bool onTrack = false;
                short TrackLen = 0;
                for (short i = 0; i < MovedPlace.size(); ++i) { // POUR CHAQUE COL IMPACTéE

                    gdetect2.clear();
                    for (short j = 0; j < ColSize; j++) {
                        if (tableau2DReplica[j][MovedPlace[i]] == INDETERMINE) {//si j arrive sur un blanc
                            if (onTrack) {//si je suis(SUIVRE) un noir
                                // FIN DE GROUPEMENT
                                onTrack = false;
                                gdetect2.push_back(TrackLen);
                                TrackLen = 0;
                            }
                        }
                        else {
                            //pas blanc donc je suis noir
                            onTrack = true;
                            TrackLen++;
                        }
                    }
                    if (TrackLen != 0) {
                        onTrack = false;
                        gdetect2.push_back(TrackLen);
                        TrackLen = 0;
                    }

                    // Version Needle 
                    signed short deltaNW = contraintesCols[MovedPlace[i]].size() - gdetect2.size();

                    std::vector<short> seq1 = contraintesCols[MovedPlace[i]];
                    std::vector<short> seq2 = gdetect2;
                    needlemanWunsch(seq1, seq2);
                    for (short j = 0; j < seq1.size(); j++) {
                        ScoresModified[i] += pow((seq1[j] - seq2[j]), 2);
                    }
                    ScoresModified[i] += pow(deltaNW, 2);

                    // OLD SCORING
                    ActualScore += ScoresCols[MovedPlace[i]];
                    TOTScoresModified += ScoresModified[i];


                }

                std::cout << std::endl << "ActualScore : " << ActualScore;
                std::cout << std::endl << "ScoresModif : " << TOTScoresModified;

                if (ActualScore >= TOTScoresModified) { // SI MIEUX, C EST TIPART
                    for (short j = 0; j < ColSize; j++) {
                        if (tableau2D[seedRow][j] != BLACK100 && tableau2D[seedRow][j] != WHITE100) {
                            tableau2D[seedRow][j] = tableau2DReplica[seedRow][j];
                        }
                    }
                    std::cout << std::endl << "Sol Accepted";
                    TrueIteration++;
                    ScoreTotal = ScoreTotal + (TOTScoresModified - ActualScore);

                    for (short i = 0; i < MovedPlace.size(); ++i) {
                        ScoresCols[MovedPlace[i]] = ScoresModified[i];
                    }


                }
                else { // SI MOINS BON
                    // Acceptance ? 
                    float rapportDerreur = ActualScore / (TOTScoresModified - ActualScore);

                    if (rapportDerreur > errorAcceptance) { // quand j ai 1 j accepte les truc aussi merqique que ma sol, si ca augmente j accepte moins 
                        for (short j = 0; j < ColSize; j++) {
                            if (tableau2D[seedRow][j] != BLACK100 && tableau2D[seedRow][j] != WHITE100) {
                                tableau2D[seedRow][j] = tableau2DReplica[seedRow][j];
                            }
                        }
                        std::cout << std::endl << "Sol Accepted";
                        TrueIteration++;
                        ScoreTotal += (TOTScoresModified - ActualScore);

                        for (short i = 0; i < MovedPlace.size(); ++i) {
                            ScoresCols[MovedPlace[i]] = ScoresModified[i];
                        }
                        errorAcceptance *= ERRORACCEPTANCEREDUCTOR;
                    }
                }
            }
        }
        else {
            // Remake Replika on real
            for (short j = 0; j < ColSize; j++) {
                tableau2DReplica[seedRow][j] = tableau2D[seedRow][j];
            }
            std::cout << std::endl << "ERROR seed want liberty";
        }
        std::cout << std::endl;
        std::cout << "scoreTotal : " << ScoreTotal;
        std::cout << std::endl;
        std::cout << "TrueIteration : " << TrueIteration;
        iteration++;
        if (ScoreTotal == 0 || iteration > 10000)RecuitRunStep = false;
        //if (ScoreTotal == 0)RecuitRunStep = false;
        //RecuitRunStep = false;
    }

    createPTOFile(cuttedNamechar, tableau2D, ScoreTotal, RowSize, ColSize);


    /*--------------------------------------------------------------------------------//Fin affichage//--------------------------------------------------------------------------------*/
    // AFFICHAGE SOL
    std::cout << std::endl << std::endl;
    std::cout << "FINAL";
    for (short i = 0; i < RowSize; i++)
    {
        std::cout << std::endl;
        for (short j = 0; j < ColSize; j++)
        {
            std::cout << tableau2D[i][j];
        }
    }

    for (short i = 0; i < RowSize; i++) {
        delete[] tableau2D[i];
    }
    delete[] tableau2D;


}//END MAIN
