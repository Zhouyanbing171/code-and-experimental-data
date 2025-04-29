import com.opencsv.CSVReader;
import com.opencsv.exceptions.CsvException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

public class Reader {
    // 读取基因本体文件，获得基因本体数据数组
    public static ArrayList<Term> readOboFile(String obo) throws IOException {
        ArrayList<Term> termList = new ArrayList<>();
        String content = new String(Files.readAllBytes(Paths.get(obo)));
        Pattern termPattern = Pattern.compile("\\[Term][\\s\\S.+?]+?\\n\\n");
        Pattern idPattern = Pattern.compile("(?<=id: )GO:\\d{7}");
        Pattern namePattern = Pattern.compile("(?<=name: ).+?\\\\n");
        Pattern namespacePattern = Pattern.compile("(?<=namespace: ).+?\\n");
        Pattern isaPattern = Pattern.compile("(?<=is_a: )GO:\\d{7}");
        Pattern partofPattern = Pattern.compile("(?<=relationship: part_of )GO:\\d{7}");
        Pattern obsoletePattern = Pattern.compile("is_obsolete: true");
        Matcher m = termPattern.matcher(content);
        while (m.find()) {
            Term termNode = new Term();
            String term = m.group();

            Matcher idm = idPattern.matcher(term);
            if (idm.find()) {
                termNode.id = idm.group().replace(":", "").trim();
            }
            Matcher names = namePattern.matcher(term);
            if (names.find()) {
                termNode.name = names.group().replace(":", "").trim();
            }
            Matcher nsm = namespacePattern.matcher(term);
            if (nsm.find()) {
                termNode.namespace = nsm.group().replace(":", "").trim();
            }
            Matcher isam = isaPattern.matcher(term);
            while (isam.find()) {
                termNode.isID.add(isam.group().replace(":", "").trim());
            }
            Matcher part_ofm = partofPattern.matcher(term);
            while (part_ofm.find()) {
                termNode.partID.add(part_ofm.group().replace(":", "").trim());
            }
            termNode.isObsolete = obsoletePattern.matcher(term).find();
            termList.add(termNode);
        }

        return termList;
    }

    // 读取注释信息文件，获得注释信息数据数组
    public static ArrayList<Annotation> readGafFile(String gaf) throws IOException {
        ArrayList<Annotation> annotationList = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(gaf))) {
            String line;
            while ((line = br.readLine()) != null) {
                // 跳过以!开头的注释行
                if (line.startsWith("!")) {
                    continue;
                }
                String[] row = line.split("\t");
                if (row.length >= 12) {
                    Annotation annotationNode = new Annotation();
                    // 第1列：标识符ID
                    annotationNode.geneID = row[2];
                    // 第2列：基因本体术语ID
                    annotationNode.goID = row[4].replace(":", "");
                    // 第3列：证据编码
                    annotationNode.evidenceCode = row[6];
                    // 第4列：基因本体术语所属分支
                    annotationNode.nameSpace = row[8];
                    // 第5列：基因本体术语的同义词
                    String synonymData = row[9];
                    annotationNode.synonym = new ArrayList<>(); // 初始化同义词列表
                    if (synonymData != null && !synonymData.isEmpty()) { // 检查是否不是空
                        String[] synonyms = synonymData.split("\\|");
                        for (String synonym : synonyms) {
                            // 过滤重复同义词
                            if (synonym != null && !synonym.trim().isEmpty()) {
                                annotationNode.synonym.add(synonym);
                            }
                        }
                    }
                    annotationList.add(annotationNode);
                }
            }
        }
        return annotationList;
    }

    // 读取基因功能网络文件，获得基因网络节点数据数组
    public static ArrayList<FunctionNet> readNetFile(String net) throws IOException {
        ArrayList<FunctionNet> netList = new ArrayList<>();
        for (String line : Files.readAllLines(Paths.get(net))) {
            FunctionNet netNode = new FunctionNet();
            netNode.gene1 = line.split("\t")[0];
            netNode.gene2 = line.split("\t")[1];
            netNode.value = Double.parseDouble(line.split("\t")[2]);

            netList.add(netNode);
        }
        return netList;
    }

    // 读取基因相似度文件，获得基因相似度数据数组 哈希表
    public static HashMap<String, Double> getSimilarityResult(String gene_type, double alpha) throws IOException {
        HashMap<String, Double> similarityResult = new HashMap<>();
        for (String line : Files.readAllLines(Paths.get("result/" + gene_type + "/similarityResult" + alpha + ".txt"))) {
            String gene1 = line.split("\t")[0];
            String gene2 = line.split("\t")[1];
            String value = line.split("\t")[2];
            similarityResult.put(gene1 + ":" + gene2, Double.parseDouble(value));
        }
        return similarityResult;
    }

    // 读取基因相似度文件，获得基因相似度数据数组 列表
    public static HashMap<String, Double> getHumanSimilarityResult(String gene_type, double alpha) throws IOException {
        // 根据映射文件获取标识符与基因的映射表
        HashMap<String, String> idGeneHashMap = new HashMap<>();
        for (String line : Files.readAllLines(Paths.get("buf/Human/id_symbol.map"))) {
            String[] values = line.split("\t");
            idGeneHashMap.put(values[1], values[0]);
        }
        // 读取基因相似度文件
        HashMap<String, Double> similarityResult = new HashMap<>();
        for (String line : Files.readAllLines(Paths.get("result/" + gene_type + "/similarityResult" + alpha + ".txt"))) {
            String gene1 = line.split("\t")[0];
            String gene2 = line.split("\t")[1];
            double value = Double.parseDouble(line.split("\t")[2]);
            similarityResult.put(idGeneHashMap.get(gene1) + ":" + idGeneHashMap.get(gene2), value);
        }
        return similarityResult;
    }

    // 读取Yeast ec文件，获得ec数据
    public static HashMap<String, HashSet<String>> getYeastEcNumGeneHashMap(ArrayList<String> filePaths) throws IOException {
        // 基因证据ec编号 哈希表 初始化
        HashMap<String, HashSet<String>> ecNumGeneHashMap = new HashMap<>();
        try (Stream<String> lines = Files.lines(Paths.get(filePaths.get(3)))) {
            lines.map(line -> line.split("\t")) // 使用制表符作为分隔符读取文件
                    .filter(values -> values.length >= 4) // 确保有足够的列
                    .forEach(values -> {
                        String ec = values[2];
                        String gene = values[3];
                        ecNumGeneHashMap.computeIfAbsent(ec, k -> new HashSet<>()).add(gene); // 自动处理缺失项
                    });
        } catch (IOException e) {
            throw new IOException(e);
        }
        // 过滤掉只有一个基因的ec编号
        ecNumGeneHashMap.entrySet().removeIf(entry -> entry.getValue().size() <= 1);
        // 移除空字符串的键（如果存在）
        ecNumGeneHashMap.remove("");
        return ecNumGeneHashMap;
    }

    // 读取Arabidopsis thaliana ec文件，获得ec数据
    public static HashMap<String, HashSet<String>> getArabidopsisEcNumGeneHashMap(ArrayList<String> filePaths) throws IOException {
        // 基因证据ec编号 哈希表 初始化
        HashMap<String, HashSet<String>> ecNumGeneHashMap = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePaths.get(3)))) {
            String line;
            // 如果有则跳过表头
            br.readLine();
            while ((line = br.readLine()) != null) {
                String[] columns = line.split("\t"); // 使用制表符作为分隔符读取文件
                if (columns.length < 8) {
                    continue;
                } // 确保行中有足够的列
                String ecNumberRaw = columns[3].trim(); // 原始EC编号
                String ecNumber = ecNumberRaw.replace("EC-", ""); // 去掉EC前缀
                String geneName = columns[6].trim(); // 基因名字
                // 过滤掉无效的基因名字
                if (geneName.equals("unknown")) {
                    continue;
                }
                // 检查EC编号和基因名字是否为空
                if (!ecNumber.isEmpty() && !geneName.isEmpty()) {
                    ecNumGeneHashMap.putIfAbsent(ecNumber, new HashSet<>()); // 如果没有则创建新的集合
                    ecNumGeneHashMap.get(ecNumber).add(geneName); // 添加基因名字到相应的EC编号
                }
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
        // 过滤掉只包含一个基因的EC编号
        ecNumGeneHashMap.values().removeIf(entry -> entry.size() <= 1);
        // 移除空字符串的键（如果存在）
        ecNumGeneHashMap.remove("");
        return ecNumGeneHashMap;
    }

    // 读取Human ec文件，获得ec数据
    public static HashMap<String, HashSet<String>> getHumanEcNumGeneHashMap(ArrayList<String> filePaths) throws IOException {
        HashMap<String, HashSet<String>> ecGeneHashMap = new HashMap<>();

        try (CSVReader csvReader = new CSVReader(new FileReader(filePaths.get(3)))) {
            List<String[]> allRows = csvReader.readAll();  // 读取所有行
            // 跳过表头
            String[] headers = allRows.get(0);

            // 找到 gene_ids 和 ecs 列的索引
            int geneIdIndex = -1;
            int ecIndex = -1;
            for (int i = 0; i < headers.length; i++) {
                if (headers[i].equalsIgnoreCase("geneids")) {
                    geneIdIndex = i;
                }
                if (headers[i].equalsIgnoreCase("ecs")) {
                    ecIndex = i;
                }
            }

            // 确保索引合法
            if (geneIdIndex == -1 || ecIndex == -1) {
                throw new IllegalArgumentException("未找到 geneids 或 ecs 列");
            }

            // 遍历余下的行
            for (int i = 1; i < allRows.size(); i++) { // 从 1 开始以跳过表头
                String[] columns = allRows.get(i);
                if (columns.length > Math.max(geneIdIndex, ecIndex)) {
                    String geneIdStr = columns[geneIdIndex].trim();
                    String ecStr = columns[ecIndex].trim();

                    if (!geneIdStr.isEmpty() && !ecStr.isEmpty()) {
                        String[] geneIdList = geneIdStr.split("\\|"); // 分割 gene_ids
                        String[] ecList = ecStr.split("\\|"); // 分割 ecs
                        // 创建对应关系
                        for (String s : ecList) {
                            String ec = s.trim();
                            // 如果EC编号为空，则跳过
                            if (ec.isEmpty()) continue;

                            // 确保基因ID能够添加到关联的EC中
                            for (String geneId : geneIdList) {
                                if (!geneId.isEmpty()) {
                                    ecGeneHashMap.computeIfAbsent(ec, k -> new HashSet<>()).add(geneId.trim()); // 添加基因 ID 到相应的 EC 编号下
                                }
                            }
                        }
                    }
                }
            }
        } catch (IOException | CsvException e) {
            throw new IOException(e); // 处理 IO 异常
        }

        // 过滤掉只包含一个EC的基因组合
        ecGeneHashMap.entrySet().removeIf(entry -> entry.getValue().size() <= 1);
        // 移除空字符串的键（如果存在）
        ecGeneHashMap.remove("");
        return ecGeneHashMap;
    }

    public static List<String> getHumanGeneNameToIdentifierMap(List<String> uniqueGenesList) throws IOException {
        List<String> geneNameToIdentifierMap = new ArrayList<>();
        // 读取 Human 基因名到 ID 的.map映射文件
        try (BufferedReader br = new BufferedReader(new FileReader("buf/Human/id_symbol.map"))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] columns = line.split("\t");
                if (columns.length > 1) {
                    String geneName = columns[0].trim();
                    String geneId = columns[1].trim();
                    if (uniqueGenesList.contains(geneName)) {
                        geneNameToIdentifierMap.add(geneId);
                    }
                }
            }
        } catch (IOException e) {
            throw new IOException(e);
        }

        return geneNameToIdentifierMap;
    }
}
